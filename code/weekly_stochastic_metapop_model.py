

"""
Spatial age-structured stochastic SEIR model with regional mixing, local mixing and built-enviroment effects.

Key Features:
- Multi-patch spatial structure with connectivity networks
- Age-structured population (2 age classes: 0-19 years, 20-65 years)
- School calendar effects on contact patterns
- Ventilation effects on transmission
- Stochastic transmission dynamics


Author: Giulia Pullano
"""

import numpy as np
import pandas as pd
from scipy.integrate import odeint
import pickle
import time
import random
import concurrent.futures
import traceback
import multiprocessing
import sys

# Set random seeds for reproducibility (None = random seed)
np.random.seed(None)
random.seed(None)


# ============================================================================
# MODEL PARAMETERS
# ============================================================================

# Fitted transmission parameter
beta = 0.153  # Base transmission rate

# Simulation settings
num_runs = 50  # Total number of simulation runs


# Population structure
age_class = 2  # Number of age classes (0: 0-19 years, 1: 20-65 years)

# Time parameters

time_steps = 28   # Number of time steps in simulation

# Compartment parameters
n_compartment = 4  # SEIR compartments (S, E, I, R)

# Disease dynamics parameters
mu = 1      # Incubation rate (E -> I transition rate)
gamma = 1   # Recovery rate (I -> R transition rate)




# ============================================================
# Load Patch Indexing
# ============================================================

codes_file = './input_data/code.csv'
data = pd.read_csv(codes_file)
n_patch = len(data)


# ============================================================
# File Paths
# ============================================================


initial_condition_file = './input_data/initial_condition_county_age_national.txt'
net_file = './input_data/network_national.txt'
school_start_file = './input_data/school_start_national.txt'
vent_file = './input_data/ventilation_weekly.csv'


# ============================================================================
# INITIALIZE DATA STRUCTURES
# ============================================================================

# Initial conditions: t0_map[patch, age_class, compartment]
t0_map = np.zeros((n_patch, age_class, n_compartment), dtype=np.int64)

# State variables over time: map_X[patch, age_class, time]
map_S = np.zeros((n_patch, age_class, time_steps + 1), dtype=np.float64)
map_E = np.zeros((n_patch, age_class, time_steps + 1), dtype=np.float64)
map_I = np.zeros((n_patch, age_class, time_steps + 1), dtype=np.float64)
map_R = np.zeros((n_patch, age_class, time_steps + 1), dtype=np.float64)

# New infections and incidence tracking
new_map_I = np.zeros((n_patch, age_class, time_steps + 1), dtype=np.float64)
incidence = np.zeros((n_patch, time_steps + 1), dtype=np.float64)
incidence2 = np.zeros((n_patch, time_steps + 1), dtype=np.float64)

# Network structure variables
degree_in = np.zeros((age_class, n_patch), dtype=np.int64)
degree_out = np.zeros((age_class, n_patch), dtype=np.int64)
prob_stay = np.zeros((age_class, n_patch))

# Effective population accounting for mobility
N_sub_eff = np.zeros(n_patch, dtype=np.float64)

# Force of infection components
lambda_v = np.zeros((age_class, n_patch))
lambda_r = np.zeros((age_class, n_patch))

# Ventilation effect data
vent = np.zeros((n_patch, time_steps + 1), dtype=np.float64)


============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

def read_initial_conditions(file):
    """
    Read initial conditions for each patch and age class.
    
    File format: patch_id, age_class, S, E, I, R (space-separated)
    
    Parameters:
    -----------
    file : str
        Path to initial conditions file
    """
    df = pd.read_csv(file, sep='\s+', header=None)
    for _, row in df.iterrows():
        idx, g_idx, S, E, I, R = row
        t0_map[int(idx)][int(g_idx)] = [S, E, I, R]


def read_vent(file):
    """
    Read ventilation effect data for each patch and time step.
    
    Ventilation values modify the transmission rate (beta).
    Missing values are forward-filled, and leading zeros are replaced
    with the first non-zero value.
    
    Parameters:
    -----------
    file : str
        Path to ventilation data file (CSV format)
    """
    df = pd.read_csv(file, sep=',')
    for _, row in df.iterrows():
        idx, var, time = row
        vent[int(idx)][int(time)] = var
    
    # Process ventilation data: fill missing values
    for idx in range(len(vent)):
        s = pd.Series(vent[idx])
        # Forward fill zeros (treat as missing)
        s = s.replace(0, np.nan).ffill()
        vent[idx] = s.fillna(0).values
        
        # Fill leading zeros with first non-zero value
        try:
            first_nonzero = vent[idx][np.nonzero(vent[idx])][0]
            first_idx = np.argmax(vent[idx] != 0)
            vent[idx][:first_idx] = first_nonzero
        except:
            continue


def read_school_start(file):
    """
    Read school start week for each patch.
    
    School calendars affect contact patterns between age groups.
    
    Parameters:
    -----------
    file : str
        Path to school start data file
        
    Returns:
    --------
    sc_start : np.array
        Array of school start weeks for each patch
    """
    df = pd.read_csv(file, sep='\s+', header=None)
    sc_start = np.zeros(n_patch, dtype=int)
    for _, row in df.iterrows():
        week_sc, idx = row
        sc_start[int(idx)] = int(week_sc)
    return sc_start


def read_connectivity_network_adult(matrix_file, num_groups):
    """
    Read adult mobility network between patches.
    
    Network represents movement probabilities between geographic patches
    for the adult age class.
    
    Parameters:
    -----------
    matrix_file : str
        Path to network file
    num_groups : int
        Number of patches
        
    Returns:
    --------
    net : np.array
        Connectivity matrix [from_patch, to_patch]
    """
    net = [[0 for i in range(num_groups)] for i in range(num_groups)]
    
    try:
        with open(matrix_file, 'r') as file:
            for line in file:
                n, m, prob = map(float, line.split())
                n, m = int(n), int(m)
                net[n][m] = prob
    except FileNotFoundError:
        print("Network file not found")
    
    return net


# ============================================================================
# LOAD DATA
# ============================================================================

# Load initial conditions
read_initial_conditions(initial_condition_file)

# Load ventilation data
read_vent(vent_file)

# Load school start dates
sc_start = read_school_start(school_start_file)

# ============================================================================
# CONTACT MATRICES
# ============================================================================

# Contact matrix when school IS in session
# c_s[age_class_i][age_class_j] = average contacts per day
c_s = {}
c_s[0] = {}
c_s[1] = {}
c_s[0][0] = 18.59  # Young-to-young contacts
c_s[0][1] = 4.21   # Young-to-adult contacts
c_s[1][0] = 5.58   # Adult-to-young contacts
c_s[1][1] = 8.84   # Adult-to-adult contacts

# Contact matrix when school is NOT in session
c_h = {}
c_h[0] = {}
c_h[1] = {}
c_h[0][0] = 7.78   # Young-to-young contacts (reduced)
c_h[0][1] = 2.55   # Young-to-adult contacts
c_h[1][0] = 5.83   # Adult-to-young contacts
c_h[1][1] = 8.15   # Adult-to-adult contacts

# Build time-varying contact matrices based on school calendars
data_array = [[{} for i in range(n_patch)] for j in range(time_steps)]
c = [[[[0 for k in range(2)] for j in range(2)] for i in range(n_patch)]
     for day in range(time_steps)]

for day in range(time_steps):
    for i in range(n_patch):
        # Use holiday contacts if before school start, else school contacts
        data_array[day][i] = c_h if day < sc_start[i] else c_s
        dd = data_array[day][i]
        c[day][i] = np.array([[dd[i][j] for j in dd[i]] for i in dd])

c = np.array(c)


# ============================================================================
# NETWORK STRUCTURE
# ============================================================================

# Load adult mobility network
num_groups = n_patch
net = read_connectivity_network_adult(net_file, num_groups)
net = np.array(net)

# Load full network to identify neighbors
tot9 = pd.read_csv('./input_data/network_national.txt', sep=' ', header=None)
nei = tot9[(tot9[2] != 0) & (tot9[1] != tot9[0])]
neig = nei.groupby(0)[1].unique().to_dict()

# Ensure every patch has an entry (even if no neighbors)
for i in range(n_patch):
    if i not in neig:
        neig[i] = []

# Create identity network for young age class (no mobility)
net_c = [[1 if i == j else 0 for j in range(num_groups)] for i in range(num_groups)]

# Combined network: net2[age_class, from_patch, to_patch]
net2 = [net_c, net]
net2 = np.array(net2)

# ============================================================================
# CALCULATE EFFECTIVE POPULATION
# ============================================================================

# Age class indices
young = 0
adult = 1

# Total population per patch and age class
t0_sum = np.sum(t0_map, axis=-1)  # Sum over SEIR compartments

# Extract population vectors
young_pop = t0_sum[:, young]  # Young population per patch
adult_pop = t0_sum[:, adult]  # Adult population per patch

# Effective population accounts for adult mobility
# N_sub_eff[i] = young_pop[i] + sum_j (net[j][i] * adult_pop[j])
N_sub_eff = young_pop + np.dot(adult_pop, net)


# ============================================================================
# SEIR SIMULATION FUNCTION
# ============================================================================

def SIER_simulation(time_steps):
    """
    Run stochastic SEIR simulation across all patches and time steps.
    
    The model includes:
    - Age-structured populations with different contact patterns
    - Spatial structure with adult mobility between patches
    - School calendar effects on contacts
    - Ventilation effects on transmission
    - Stochastic transitions between compartments
    
    Parameters:
    -----------
    time_steps : int
        Number of time steps to simulate
        
    Returns:
    --------
    inc : np.array
        Incidence (new infections) per patch and time step
    """
    
    # Initialize compartments from initial conditions
    S = t0_map[:, :, 0].copy()
    E = t0_map[:, :, 1].copy()
    I = t0_map[:, :, 2].copy()
    R = t0_map[:, :, 3].copy()
    
    # Initialize tracking matrices
    map_S = np.zeros((n_patch, age_class, time_steps))
    map_E = np.zeros((n_patch, age_class, time_steps))
    map_I = np.zeros((n_patch, age_class, time_steps))
    map_R = np.zeros((n_patch, age_class, time_steps))
    inc = np.zeros((n_patch, time_steps))
    lambda_ = np.zeros((n_patch, time_steps))
    
    # Set initial values
    map_S[:, :, 0] = S
    map_E[:, :, 0] = E
    map_I[:, :, 0] = I
    map_R[:, :, 0] = R
    inc[:, 0] = np.sum(E, axis=1)
    
    # Main simulation loop
    for day in range(time_steps):
        I_day = map_I[:, :, day]
        
        for i in range(n_patch):
            N_eff_i = N_sub_eff[i]
            h_arr = np.array(neig[i])  # Neighbor patches
            
            for g in range(age_class):
                # Calculate force of infection (FOI)
                # Local transmission within patch i
                foi = (net2[g, i, i] *
                       np.sum(net2[:, :, i] * I_day.T *
                              c[day, i, g, :].reshape(-1, 1), axis=(0, 1)) / N_eff_i)
                
                # Add mobility-related transmission for adults
                if (g == 1) and len(h_arr) > 0:
                    # Infections from adults visiting from neighbor patches
                    I_sub_eff = np.sum(net2[g, :n_patch, h_arr] *
                                      map_I[:n_patch, 1, day], axis=1)
                    foi += np.sum(net2[g, i, h_arr] *
                                 c[day, h_arr, g, 1] *
                                 I_sub_eff / N_sub_eff[h_arr])
                
                # Apply ventilation effect to transmission rate
                beta2 = beta + vent[i][day]
                if beta2 < 0:
                    beta2 = 0
                
                lambda_[i, day] = beta2 * foi
                
                # Get current compartment values
                S_ig = map_S[i, g, day]
                E_ig = map_E[i, g, day]
                I_ig = map_I[i, g, day]
                R_ig = map_R[i, g, day]
                
                # Stochastic transitions between compartments
                # S -> E: New exposures
                prob_infection = 1 - np.exp(-lambda_[i, day])
                new_exposed = (np.random.binomial(int(S_ig), prob_infection)
                              if S_ig > 0 else 0)
                
                # E -> I: Progression to infectious
                new_infections = (np.random.binomial(int(E_ig), mu)
                                 if E_ig > 0 else 0)
                
                # I -> R: Recovery
                new_recovered = (np.random.binomial(int(I_ig), gamma)
                                if I_ig > 0 else 0)
                
                # Update compartments for next time step
                map_S[i, g, day + 1] = S_ig - new_exposed
                map_E[i, g, day + 1] = E_ig - new_infections + new_exposed
                map_I[i, g, day + 1] = I_ig + new_infections - new_recovered
                map_R[i, g, day + 1] = R_ig + new_recovered
                
                # Track incidence (new infections)
                inc[i, day + 1] += new_infections
    
    return inc


# ============================================================================
# RUN SIMULATIONS AND SAVE RESULTS
# ============================================================================

print(f"Starting {num_runs} simulation runs...")
print(f"Parameters: beta={beta}, gamma={gamma}, mu={mu}")

for run in range(num_runs):
    print(f"Running simulation {run + 1}/{num_runs}...", end='\r')
    
    # Run simulation
    inc = SIER_simulation(time_steps)
    
    # Save results
    output_file = (f'./output/fit/simulated_incidence_run_{run}_'
                   f'beta_{str(beta)[:5]}_gamma_{gamma}_mu_{mu}_v4.csv')
    pd.DataFrame(inc).T.to_csv(output_file)

print(f"\nSimulations complete! Results saved to ./output/")
