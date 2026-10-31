# 05 — Materials Characterization Workbench

Data-processing utilities for illustrative XRD, FTIR, DSC, and TGA signals.

- XRD: local peak candidates and d-spacing from Bragg's law.
- FTIR: local maxima in an absorbance spectrum.
- DSC: temperature of the largest deviation from the median heat-flow baseline.
- TGA: percentage mass loss.

Run: `python characterization.py`

Peak detection is sensitive to background, noise, resolution, and thresholds. Candidate peaks do not establish phase identity or chemical assignment. DSC event sign depends on instrument convention. Example signals are synthetic.
