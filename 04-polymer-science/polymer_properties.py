"""Polymer molecular-weight averages and Fox blend Tg estimate."""
import numpy as np


def molecular_weight_averages(molecular_weights, molecule_counts):
    m, n = np.asarray(molecular_weights, float), np.asarray(molecule_counts, float)
    if m.ndim != 1 or n.ndim != 1 or m.size == 0 or m.size != n.size:
        raise ValueError("Equal-length non-empty 1D arrays are required.")
    if np.any(m <= 0) or np.any(n < 0) or n.sum() <= 0:
        raise ValueError("Molecular weights must be positive and counts non-negative with positive total.")
    mn = np.sum(n * m) / np.sum(n)
    mw = np.sum(n * m * m) / np.sum(n * m)
    return {"Mn": float(mn), "Mw": float(mw), "dispersity": float(mw / mn)}


def fox_blend_tg(weight_fractions, component_tg_k):
    w, tg = np.asarray(weight_fractions, float), np.asarray(component_tg_k, float)
    if w.ndim != 1 or tg.ndim != 1 or not w.size or w.size != tg.size:
        raise ValueError("Equal-length non-empty 1D arrays are required.")
    if np.any(w < 0) or np.any(tg <= 0) or not np.isclose(w.sum(), 1.0, atol=1e-8):
        raise ValueError("Fractions must be non-negative and sum to 1; Tg must be in kelvin.")
    return float(1.0 / np.sum(w / tg))


if __name__ == "__main__":
    print(molecular_weight_averages([10_000, 30_000, 80_000], [50, 35, 15]))
    tg = fox_blend_tg([0.6, 0.4], [350, 410])
    print(f"Blend Tg: {tg:.2f} K ({tg - 273.15:.2f} °C)")
