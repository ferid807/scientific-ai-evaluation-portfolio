# 03 — Analytical Chemistry: Calibration & Quality

Fits a linear calibration curve and calculates residuals, R², and illustrative LOD/LOQ using the explicitly chosen convention LOD = 3.3 σblank / slope and LOQ = 10 σblank / slope.

Run: `python calibration_analysis.py`

R² alone does not establish calibration validity. Inspect residuals, range, heteroscedasticity, weighting, replicates, matrix effects, selectivity, and method-specific acceptance criteria. LOD/LOQ conventions vary by method and standard. Example responses are synthetic.
