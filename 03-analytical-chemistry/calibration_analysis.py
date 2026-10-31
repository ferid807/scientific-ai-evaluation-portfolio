"""Linear calibration, residuals and an explicitly stated LOD/LOQ convention."""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class CalibrationResult:
    slope: float
    intercept: float
    r_squared: float
    residuals: np.ndarray
    lod: float
    loq: float


def fit_calibration(concentration, response, blank_sd):
    x, y = np.asarray(concentration, float), np.asarray(response, float)
    if x.ndim != 1 or y.ndim != 1 or x.size != y.size or x.size < 3:
        raise ValueError("Equal-length 1D arrays with at least three points are required.")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or blank_sd < 0 or np.ptp(x) == 0:
        raise ValueError("Inputs must be finite, blank SD non-negative, and x must span a range.")
    slope, intercept = np.polyfit(x, y, 1)
    if slope <= 0:
        raise ValueError("Expected a positive calibration slope.")
    residuals = y - (slope * x + intercept)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot else float("nan")
    return CalibrationResult(float(slope), float(intercept), float(r2), residuals,
                             float(3.3 * blank_sd / slope), float(10.0 * blank_sd / slope))


def estimate_concentration(response, slope, intercept, response_sd=0.0):
    if slope <= 0 or response_sd < 0:
        raise ValueError("Slope must be positive and response SD non-negative.")
    return (response - intercept) / slope, response_sd / slope


if __name__ == "__main__":
    # Synthetic illustrative instrument responses.
    result = fit_calibration([0, 1, 2, 3, 4, 5], [0.04, 0.51, 1.03, 1.48, 2.02, 2.49], 0.025)
    print(f"Response = {result.slope:.4f} C + {result.intercept:.4f}")
    print(f"R²={result.r_squared:.5f}; LOD={result.lod:.4f}; LOQ={result.loq:.4f}")
