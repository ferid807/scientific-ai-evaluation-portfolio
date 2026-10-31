# 04 — Polymer Science Toolkit

Calculates number-average molecular weight (Mn), weight-average molecular weight (Mw), dispersity (Đ), and an idealized blend Tg using the Fox equation.

```
Mn = sum(Ni Mi) / sum(Ni)
Mw = sum(Ni Mi^2) / sum(Ni Mi)
Đ = Mw / Mn
1/Tg(blend) = sum(wi/Tg,i)
```

Run: `python polymer_properties.py`

Fox equation inputs use kelvin and mass fractions summing to one. It is an idealized estimate for miscible blends; immiscibility, phase separation, specific interactions, and thermal history can cause deviations.
