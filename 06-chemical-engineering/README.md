# 06 — Chemical Engineering Calculation Suite

Basic steady-state mass balance, sensible heat duty, and ideal first-order PFR/CSTR sizing.

```
mass residual = total inlet - total outlet
Q = m_dot Cp deltaT
V_PFR = (q/k) ln[1/(1-X)]
V_CSTR = q X / [k(1-X)]
```

Run: `python calculations.py`

Reactor formulas assume steady state, constant density/volumetric flow, isothermal operation, ideal reactor behavior, and irreversible first-order kinetics. Actual design also requires hydraulics, pressure drop, heat transfer, safety, controls, phase behavior, and validated physical properties.
