# Scientific AI Evaluation & Scientific Computing Portfolio

Reproducible Python projects spanning thermodynamics, reaction kinetics, analytical chemistry, polymer science, materials characterization, and chemical engineering calculations.

This portfolio demonstrates scientific computing and the use of independently checkable calculations, tests, assumptions, and evidence to evaluate AI-generated technical explanations.

> **Data integrity:** Example signals and datasets are synthetic/illustrative. They are not experimental measurements or results from an industrial laboratory.

## Projects

| Project | Scope | Main methods |
|---|---|---|
| [01 — Thermodynamics & Equilibrium](01-thermodynamics-equilibrium/) | Ideal-gas equilibrium for A + B ⇌ C | Reaction quotient, extent of reaction, van't Hoff relation |
| [02 — Reaction Kinetics](02-reaction-kinetics/) | First-order batch reaction profiles | Arrhenius equation, analytical and numerical ODE solutions |
| [03 — Analytical Chemistry](03-analytical-chemistry/) | Calibration and measurement quality | Regression, residuals, illustrative LOD/LOQ |
| [04 — Polymer Science](04-polymer-science/) | Molecular-weight distributions and blend Tg | Mn, Mw, dispersity, Fox equation |
| [05 — Materials Characterization](05-materials-characterization/) | XRD, FTIR, DSC and TGA signal utilities | Bragg's law, peak detection, thermal metrics |
| [06 — Chemical Engineering](06-chemical-engineering/) | Balances and ideal reactor sizing | Conservation equations, PFR/CSTR design |

## Quick start

Python 3.10+ is recommended.

```bash
git clone https://github.com/ferid807/scientific-ai-evaluation-portfolio.git
cd scientific-ai-evaluation-portfolio
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
python -m pytest
```

Run an individual example from its folder, e.g. `python 02-reaction-kinetics/kinetics.py`.

## Structure

```text
01-thermodynamics-equilibrium/
02-reaction-kinetics/
03-analytical-chemistry/
04-polymer-science/
05-materials-characterization/
06-chemical-engineering/
scientific-ai-evaluation/
tests/
```

See [scientific-ai-evaluation/evaluation_rubric.md](scientific-ai-evaluation/evaluation_rubric.md) for a practical review rubric. Example evaluations are illustrative, not an externally validated benchmark.

## Engineering principles

- Prefer SI units internally and label units in all inputs/outputs.
- State assumptions and applicability limits.
- Compare numerical solutions with analytical results where possible.
- Make synthetic data unmistakably synthetic.
- Do not claim experimental validation without suitable experimental evidence.

These educational tools are not certified engineering design software, laboratory protocols, or substitutes for professional review. Real applications may require non-ideal thermodynamics, validated property data, transport phenomena, detailed mechanisms, safety analysis, and experimental validation.

## Author

**Farid Karimli** — Chemical Engineering | Materials & Polymer Science | Scientific Computing

[GitHub profile](https://github.com/ferid807)
