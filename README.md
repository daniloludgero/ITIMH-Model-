# ITIMH-Model
# ITIMH Framework – Figure Generation Code

This repository contains Python code to reproduce **Figures 2–7** of the manuscript:

> **Induced Trophic-Immunological Mimicry Hypothesis (ITIMH): A Unified Dynamic Immunomechanical Framework with Empirical Validation Against Published Biomechanical, Transcriptomic, and Clinical Data**  
> *Danilo Ludgero (Independent Researcher)*

The code implements the mathematical models described in the paper:  
- Closed‑form elimination probability  
- Dynamic damage accumulation with tumour recovery (`γ`)  
- ODE kinetics of YAP/TAZ nuclear translocation  
- Absorbing Markov chain population dynamics  
- Bayesian variance decomposition (visualisation only)  

All parameters are taken directly from the manuscript and its cited literature (Cross 2007, Meng 2016, Zhao 2010, etc.).

---

## Requirements

Install the following Python packages:

```bash
pip install numpy matplotlib scipy