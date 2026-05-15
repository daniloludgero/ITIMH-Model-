
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
```

The code was tested with:

· Python 3.8+
· NumPy 1.21+
· Matplotlib 3.5+
· SciPy 1.7+

---

File contents

File Description
generate_figures_ITIMH.py Main script that generates all figures (2 to 7)
README.md This file

---

How to run

1. Save the script as generate_figures_ITIMH.py
2. Open a terminal in the same folder
3. Run:

```bash
python generate_figures_ITIMH.py
```

The script will:

· Print Pearson correlation coefficients for Fig. 2 and Fig. 6 (as a sanity check)
· Display each figure interactively (pop‑up windows)
· Save all figures as PNG (300 dpi) and PDF (vector) in the current directory:
  · Figura2_phase_ORR.png / .pdf
  · Figura3_AFM_lambdaTc.png / .pdf
  · Figura4_YAP_TAZ_ODE.png / .pdf
  · Figura5_recovery_kinetics.png / .pdf
  · Figura6_FAK_trials.png / .pdf
  · Figura7_Markov_dynamics.png / .pdf

---

What each figure reproduces (with exact data from the manuscript)

Figure Description Data source in paper
Fig. 2 Phase diagram (P_elim vs λTc) and correlation with published clinical ORR (r=0.985) Test 1 (p.7): λTc and ORR for 7 tumour types
Fig. 3 AFM Young’s moduli → λTc; cancer vs normal cells; power‑law E→λ Test 2 (p.7): values 7.21±4.39 (cancer) and 0.71±0.10 (normal)
Fig. 4 YAP/TAZ ODE kinetics (stiff/soft substrates) compared to Meng 2016 live‑cell data; mapping to λ and predicted drug effect Test 3 (p.8–9): rate constants from Zhao 2010, Meng 2016
Fig. 5 Dynamic damage accumulation: sustained vs fragmented contacts; killing probability vs recovery rate γ Test 4 (p.9–10): γ = 0.15 min⁻¹, 50% reduction for fragmented contacts
Fig. 6 FAK inhibitor clinical trials (ORR baseline vs combination); correlation with model predictions (r=0.829); λ variance decomposition Test 5 (p.10–11): NSCLC, pancreatic, ovarian data
Fig. 7 Absorbing Markov chain (S0→S1→S2); elimination kinetics; time to 50% elimination; S1 peak occupancy Test 6 (p.11): transition probabilities for baseline, ITIMH, ITIMH+λ‑reduction

---

Notes on reproducibility

· The script does not require any external data files – all numerical values are hard‑coded from the manuscript.
· Experimental YAP/TAZ time points (Test 3) were digitised approximately from Meng et al. 2016, Figure 2. Small differences between generated and published plots are expected but the trends and RMSE (41.7%) match the paper.
· For the power‑law conversion in Test 2 (E → λ), the exponent β = 0.92 and reference modulus E0 = 10 kPa are used as stated in Table 2 (p.5).
· The FAK trial correlation (Fig. 6) reproduces the reported r = 0.829 using the ORR values from the manuscript.

---

License & reuse

This code is provided under CC‑BY 4.0 (same as the preprint). You are free to use, adapt, and redistribute it with appropriate credit to the original article.

---

Contact

Danilo Ludgero – independent researcher
Correspondence: danilo‑ben10@hotmail.com
Preprint DOI: https://doi.org/10.5281/zenodo.20213912
