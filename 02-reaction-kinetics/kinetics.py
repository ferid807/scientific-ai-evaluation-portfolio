"""First-order batch kinetics and Arrhenius calculations."""
import numpy as np
from scipy.integrate import solve_ivp

R = 8.314462618


def arrhenius_rate_constant(pre_exponential_s, activation_energy_j_mol, temperature_k):
    if pre_exponential_s <= 0 or temperature_k <= 0 or activation_energy_j_mol < 0:
        raise ValueError("A and T must be positive; Ea cannot be negative.")
    return float(pre_exponential_s * np.exp(-activation_energy_j_mol / (R * temperature_k)))


def first_order_profile(c0_mol_l, k_s, time_s):
    t = np.asarray(time_s, dtype=float)
    if c0_mol_l < 0 or k_s < 0 or np.any(t < 0):
        raise ValueError("Concentration, k, and time must be non-negative.")
    return c0_mol_l * np.exp(-k_s * t), 1.0 - np.exp(-k_s * t)


def integrate_batch_first_order(c0_mol_l, k_s, time_s):
    t = np.asarray(time_s, dtype=float)
    if c0_mol_l < 0 or k_s < 0 or t.ndim != 1 or not len(t) or np.any(t < 0) or np.any(np.diff(t) < 0):
        raise ValueError("Supply non-negative inputs and a sorted non-empty 1D time array.")
    if t[-1] == 0:
        return np.full_like(t, c0_mol_l)
    solution = solve_ivp(lambda _, y: [-k_s * y[0]], (0.0, float(t[-1])), [c0_mol_l],
                         t_eval=t, rtol=1e-9, atol=1e-12)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y[0]


if __name__ == "__main__":
    t = np.linspace(0, 300, 101)
    ca, conversion = first_order_profile(1.0, 0.01, t)
    print(f"At 300 s: CA={ca[-1]:.4f} mol/L; conversion={conversion[-1]:.2%}")
