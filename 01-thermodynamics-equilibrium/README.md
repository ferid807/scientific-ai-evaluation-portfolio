# 01 — Thermodynamics & Chemical Equilibrium

Ideal-gas equilibrium for the model reaction A + B ⇌ C. The solver calculates the reaction quotient and equilibrium composition using reaction extent. A van't Hoff utility estimates temperature dependence assuming constant reaction enthalpy.

For ideal gases with dimensionless activities aᵢ = pᵢ/p°:

```
Qp = (pC/p°) / [(pA/p°)(pB/p°)]
ln(K2/K1) = -(ΔH°/R) × (1/T2 - 1/T1)
```

Run: `python equilibrium_solver.py`

Assumptions: one reaction, ideal-gas activities, 1:1:1 stoichiometry, and constant ΔH° in the van't Hoff approximation. This is not a general chemical-equilibrium package.
