# Human Behavior and Built Environments Drive Influenza Seasonality in the US

## Authors

- **Giulia Pullano** – Department of Biology, Georgetown University, Washington, DC, USA  
- **Andrew Tiu** – Department of Biology, Georgetown University, Washington, DC, USA  
- **Jelena Srebric** – Department of Mechanical Engineering, University of Maryland, College Park, MD, USA  
- **Donald K. Milton** – School of Public Health, University of Maryland, College Park, MD, USA  
- **Linsey C. Marr** – Department of Civil and Environmental Engineering, Virginia Tech, Blacksburg, VA, USA  
- **Shweta Bansal** – Department of Biology, Georgetown University, Washington, DC, USA  

---

## Overview

This repository contains the code (and processed datasets, to be added) supporting the manuscript:

**Human Behavior and Built Environments Drive Influenza Seasonality in the US.**

The study investigates how human behavioral patterns and built environment characteristics jointly shape the timing and spatial structure of influenza seasonality across the United States.

Specifically, we integrate:

- School calendar timing  
- Human mobility networks  
- Indoor ventilation dynamics  
- Weather variability  
- Housing characteristics  
- Social vulnerability indicators  

These drivers are incorporated into a spatially explicit, age-structured stochastic SEIR metapopulation model to quantify their relative contributions to influenza seasonal dynamics.

The framework enables evaluation of how behavioral and environmental mechanisms modulate epidemic timing, synchrony, and geographic diffusion.

---

## Code Availability

### Transmission Model


We developed a stochastic, spatially explicit, age-structured SEIR metapopulation model at the U.S. county level to mechanistically evaluate drivers of influenza seasonality. The model incorporates three transmission pathways—local transmission within counties, importation from visiting infectious individuals, and infection acquired during travel and reintroduced upon return—using county-level mobility networks and age-specific contact matrices (two age classes: 5–18 and 18–65). School-term forcing modulates age-specific contact rates, and the baseline transmission parameter (β₀) was calibrated to reproduce observed national incidence patterns (weeks 31–6). 

The model was initialized using observed county-level summer incidence to preserve realistic spatial heterogeneity. We systematically evaluated mechanistic contributions by sequentially adding (i) spatial connectivity, (ii) connectivity with school-term forcing, and (iii) local transmissibility modulation capturing crowding and ventilation effects through a residual transmission enhancement framework. Model performance was assessed based on its ability to reproduce observed epidemic onset timing patterns (2016/17–2019/20 seasons).


### Multi-regression Model

We implemented generalized linear regression models to quantify the associations between summer influenza incidence and epidemic onset timing with hypothesized seasonal drivers, including regional mixing, indoor social mixing, and built-environment characteristics. All predictors were z-normalized, multicollinearity was assessed using variance inflation factors (VIF > 5 excluded), and Driscoll–Kraay heteroskedasticity- and autocorrelation-consistent standard errors were used to account for spatial and temporal dependence across U.S. counties (2016/17–2019/20 seasons).


## Citation

If you use this code, please cite:

Pullano G, Tiu A, Srebric J, Milton DK, Marr LC, Bansal S.  
*Human Behavior and Built Environments Drive Influenza Seasonality in the US.*

(Full citation will be updated upon publication.)
