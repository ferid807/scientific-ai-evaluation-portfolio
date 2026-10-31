"""Basic peak detection and thermal characterization utilities."""
import numpy as np
from scipy.signal import find_peaks


def bragg_d_spacing(two_theta_deg, wavelength, order=1):
    if wavelength <= 0 or order <= 0 or not 0 < two_theta_deg < 180:
        raise ValueError("Require wavelength/order > 0 and 0 < 2theta < 180 degrees.")
    return float(order * wavelength / (2 * np.sin(np.deg2rad(two_theta_deg / 2))))


def xrd_peak_positions(two_theta_deg, intensity, prominence=None):
    x, y = np.asarray(two_theta_deg, float), np.asarray(intensity, float)
    if x.ndim != 1 or y.ndim != 1 or x.size != y.size or x.size < 3 or np.any(np.diff(x) <= 0):
        raise ValueError("Equal-length 1D arrays and strictly increasing x are required.")
    peaks, properties = find_peaks(y, prominence=prominence)
    return x[peaks], y[peaks], properties


def ftir_absorption_bands(wavenumber_cm1, absorbance, prominence=0.01):
    return xrd_peak_positions(wavenumber_cm1, absorbance, prominence)[:2]


def mass_loss_percent(initial_mass, final_mass):
    if initial_mass <= 0 or final_mass < 0 or final_mass > initial_mass:
        raise ValueError("Require initial mass > 0 and 0 <= final mass <= initial mass.")
    return 100.0 * (initial_mass - final_mass) / initial_mass


def dsc_extremum_temperature(temperature_c, heat_flow):
    t, h = np.asarray(temperature_c, float), np.asarray(heat_flow, float)
    if t.ndim != 1 or h.ndim != 1 or t.size != h.size or t.size < 3 or np.any(np.diff(t) <= 0):
        raise ValueError("Equal-length 1D arrays and increasing temperature are required.")
    return float(t[np.argmax(np.abs(h - np.median(h)))])


if __name__ == "__main__":
    two_theta = np.linspace(10, 80, 701)
    intensity = 5 + 100 * np.exp(-0.5 * ((two_theta - 30) / 0.35) ** 2)
    peaks, _, _ = xrd_peak_positions(two_theta, intensity, prominence=10)
    print("Synthetic XRD peak 2theta:", peaks)
    print("d-spacing (Cu Kalpha, angstrom):", bragg_d_spacing(float(peaks[0]), 1.5406))
    print("Illustrative TGA mass loss (%):", mass_loss_percent(10, 7.2))
