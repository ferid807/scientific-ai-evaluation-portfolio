import importlib.util
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_equilibrium_quotient_closes():
    m = load("equilibrium", Path("01-thermodynamics-equilibrium/equilibrium_solver.py"))
    result = m.equilibrium_from_extent({"A": 1, "B": 1, "C": 0}, 1, 4)
    assert result.reaction_quotient == pytest.approx(4, rel=1e-7)
    assert min(result.moles.values()) >= 0


def test_kinetics_numerical_matches_analytical():
    m = load("kinetics", Path("02-reaction-kinetics/kinetics.py"))
    t = np.linspace(0, 100, 11)
    analytical, _ = m.first_order_profile(2, 0.03, t)
    numerical = m.integrate_batch_first_order(2, 0.03, t)
    assert np.allclose(analytical, numerical, rtol=1e-6, atol=1e-8)


def test_calibration_regression():
    m = load("calibration", Path("03-analytical-chemistry/calibration_analysis.py"))
    result = m.fit_calibration([0, 1, 2, 3], [1, 3, 5, 7], 0.1)
    assert result.slope == pytest.approx(2)
    assert result.intercept == pytest.approx(1)
    assert result.r_squared == pytest.approx(1)


def test_polymer_averages_and_fox():
    m = load("polymer", Path("04-polymer-science/polymer_properties.py"))
    result = m.molecular_weight_averages([1000, 2000], [1, 1])
    assert result["Mn"] == pytest.approx(1500)
    assert result["Mw"] == pytest.approx(1666.6666667)
    assert m.fox_blend_tg([0.5, 0.5], [300, 400]) == pytest.approx(342.8571429)


def test_materials_characterization():
    m = load("characterization", Path("05-materials-characterization/characterization.py"))
    assert m.bragg_d_spacing(30, 1.5) == pytest.approx(1.5 / (2 * np.sin(np.deg2rad(15))))
    assert m.mass_loss_percent(10, 7) == pytest.approx(30)


def test_chemical_engineering_calculations():
    m = load("engineering", Path("06-chemical-engineering/calculations.py"))
    assert m.mass_balance_residual([2, 3], [4.8]) == pytest.approx(0.2)
    assert m.sensible_heat_duty(1, 1000, 10) == pytest.approx(10000)
    assert m.pfr_volume_first_order(0.01, 0.2, 0.8) < m.cstr_volume_first_order(0.01, 0.2, 0.8)
