"""Ideal-gas equilibrium for A + B <-> C, solved by bisection on extent."""
from dataclasses import dataclass
import math

R = 8.314462618  # J mol-1 K-1

@dataclass(frozen=True)
class EquilibriumResult:
    extent_mol: float
    moles: dict
    partial_pressures_bar: dict
    reaction_quotient: float
    equilibrium_constant: float


def reaction_quotient(partial_pressures_bar, standard_pressure_bar=1.0):
    """Dimensionless Qp for ideal gases; pressures and standard pressure use bar."""
    if standard_pressure_bar <= 0:
        raise ValueError("Standard pressure must be positive.")
    pa, pb, pc = (float(partial_pressures_bar.get(s, 0.0)) for s in ("A", "B", "C"))
    if min(pa, pb, pc) < 0:
        raise ValueError("Partial pressures cannot be negative.")
    if pa == 0 or pb == 0:
        return math.inf if pc > 0 else 0.0
    return (pc / standard_pressure_bar) / ((pa / standard_pressure_bar) * (pb / standard_pressure_bar))


def vanthoff_k(k_ref, t_ref_k, t_k, delta_h_j_mol):
    """K(T) from constant-enthalpy van't Hoff approximation."""
    if k_ref <= 0 or t_ref_k <= 0 or t_k <= 0:
        raise ValueError("K and temperatures must be positive.")
    return k_ref * math.exp(-delta_h_j_mol / R * (1.0 / t_k - 1.0 / t_ref_k))


def equilibrium_from_extent(initial_moles, total_pressure_bar, k_eq, standard_pressure_bar=1.0, tolerance=1e-10):
    if total_pressure_bar <= 0 or k_eq <= 0 or standard_pressure_bar <= 0:
        raise ValueError("Pressures and equilibrium constant must be positive.")
    a0, b0, c0 = (float(initial_moles.get(s, 0.0)) for s in ("A", "B", "C"))
    if min(a0, b0, c0) < 0:
        raise ValueError("Initial moles cannot be negative.")
    lo, hi = -c0, min(a0, b0)

    def state(x):
        n = {"A": a0 - x, "B": b0 - x, "C": c0 + x}
        nt = sum(n.values())
        if nt <= 0:
            raise ValueError("Total moles must remain positive.")
        p = {s: total_pressure_bar * amount / nt for s, amount in n.items()}
        return n, p, reaction_quotient(p, standard_pressure_bar)

    def residual(x):
        return state(x)[2] - k_eq

    flo, fhi = residual(lo), residual(hi)
    if flo == 0:
        xi = lo
    elif fhi == 0:
        xi = hi
    elif flo * fhi > 0:
        raise ValueError("No equilibrium root in the feasible extent interval.")
    else:
        for _ in range(300):
            mid = (lo + hi) / 2
            fm = residual(mid)
            if abs(fm) < tolerance or abs(hi - lo) < tolerance:
                lo = hi = mid
                break
            if flo * fm <= 0:
                hi, fhi = mid, fm
            else:
                lo, flo = mid, fm
        xi = (lo + hi) / 2
    n, p, q = state(xi)
    return EquilibriumResult(xi, n, p, q, k_eq)


if __name__ == "__main__":
    result = equilibrium_from_extent({"A": 1.0, "B": 1.0, "C": 0.0}, 1.0, 4.0)
    print(f"Extent: {result.extent_mol:.6f} mol")
    print("Equilibrium moles:", result.moles)
    print(f"Qp={result.reaction_quotient:.6f}; Kp={result.equilibrium_constant:.6f}")
