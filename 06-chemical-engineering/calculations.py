"""Basic conservation and ideal first-order reactor sizing calculations."""
import math


def mass_balance_residual(inlet_mass_flows, outlet_mass_flows):
    inlet, outlet = list(map(float, inlet_mass_flows)), list(map(float, outlet_mass_flows))
    if not inlet or not outlet or min(inlet + outlet) < 0:
        raise ValueError("Provide non-empty lists of non-negative mass flows.")
    return sum(inlet) - sum(outlet)


def sensible_heat_duty(mass_flow_kg_s, heat_capacity_j_kg_k, delta_t_k):
    if mass_flow_kg_s < 0 or heat_capacity_j_kg_k <= 0:
        raise ValueError("Mass flow must be non-negative and heat capacity positive.")
    return mass_flow_kg_s * heat_capacity_j_kg_k * delta_t_k


def _validate(flow, k, conversion):
    if flow <= 0 or k <= 0 or not 0 <= conversion < 1:
        raise ValueError("Flow and k must be positive; conversion must be in [0, 1).")


def pfr_volume_first_order(flow_m3_s, k_s, conversion):
    _validate(flow_m3_s, k_s, conversion)
    return flow_m3_s / k_s * math.log(1.0 / (1.0 - conversion))


def cstr_volume_first_order(flow_m3_s, k_s, conversion):
    _validate(flow_m3_s, k_s, conversion)
    return flow_m3_s * conversion / (k_s * (1.0 - conversion))


if __name__ == "__main__":
    print("Mass balance residual (kg/s):", mass_balance_residual([2, 3], [4.8]))
    print("Sensible heat duty (W):", sensible_heat_duty(1.2, 4180, 25))
    print("PFR volume (m3):", pfr_volume_first_order(0.01, 0.2, 0.8))
    print("CSTR volume (m3):", cstr_volume_first_order(0.01, 0.2, 0.8))
