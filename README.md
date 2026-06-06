# sptPALM analysis of CAZ proteins

<p align="center">
  <img src="figures/readme/trajectory.gif" width="850">
</p>

Computational pipeline for single-particle tracking PALM (sptPALM) analysis of synaptic protein mobility and membrane dynamics under different cellular perturbations.

---

## Welcome!

Single-particle tracking photoactivated localization microscopy (sptPALM) enables quantitative analysis of membrane protein dynamics at the single-molecule level. This repository contains a reproducible computational pipeline developed for processing and analysis of sptPALM trajectories of cytomatrix active zone (CAZ) proteins under oxidative stress conditions. The workflow includes trajectory preprocessing, diffusion analysis, velocity-based mobility characterization, Hidden Markov Model inference, sensitivity analysis, and automated figure generation

The repository is actively evolving together with the project itself. New analysis modules, visualization approaches, statistical methods, and trajectory-processing tools will continue to be added over time.

---

## Project overview

The pipeline was developed for quantitative analysis of SNAP25 and Syntaxin mobility in HEK293T cells under control and oxidative stress conditions.

Experimental conditions included:
- Control
- 500 μM H₂O₂
- 500 μM H₂O₂ (30 min exposure)

The workflow processes TrackMate trajectory exports and produces diffusion, velocity, mobility-state, and robustness analyses together with publication-ready figures.
---

## Biological background

Synaptic membrane proteins exhibit heterogeneous diffusion behavior due to membrane organization, transient confinement, and protein interactions.

Single-particle tracking PALM (sptPALM) enables nanoscale characterization of these dynamics by reconstructing trajectories of individual molecules. Reliable biological interpretation, however, requires careful trajectory filtering, artifact detection, and reproducible computational analysis.
---
## Reproducibility

The analysis workflow is organized as a sequence of independent notebooks covering preprocessing, diffusion analysis, velocity analysis, mobility-state inference, sensitivity testing, and final figure generation.

All intermediate results are exported as CSV tables and can be inspected independently. The final figures are generated automatically from processed results, enabling full reproducibility of the analytical workflow.

---

## Experimental setup

### Biological system
- HEK293T cells
- SNAP25 and Syntaxin expression

### Conditions
- Control
- 500 μM H₂O₂
- 500 μM H₂O₂, 30 min exposure

### Imaging
- TIRF microscopy
- sptPALM acquisition
- 20 ms exposure
- 405 / 488 / 561 nm lasers

### Initial trajectory reconstruction
- TrackMate (Fiji)

---

## Pipeline overview

### `01_trackmate_preprocessing.ipynb`

Initial preprocessing of TrackMate trajectory exports.

This notebook performs:
- loading of TrackMate CSV trajectory files,
- metadata extraction,
- trajectory quality control,
- edge-artifact filtering,
- minimum track-length filtering,
- preprocessing statistics generation,
- export of cleaned trajectory tables.

Generated outputs include:
- cleaned localization tables,
- preprocessing QC summaries,
- edge-filter diagnostics.

---

### `02_diffusion_analysis.ipynb`

Mean squared displacement (MSD) and diffusion analysis.

This notebook performs:
- MSD calculation,
- estimation of diffusion coefficients using

```text
D = slope / 4
```

for two-dimensional diffusion,
- per-track diffusion analysis,
- per-cell aggregation,
- diffusion distribution visualization.

Generated outputs include:
- MSD tables,
- diffusion coefficient summaries,
- diffusion figures.

---

### `03_velocity_analysis.ipynb`

Frame-to-frame displacement and velocity analysis.

This notebook performs:
- frame-to-frame displacement calculation,
- velocity estimation,
- velocity distribution analysis,
- empirical and Rayleigh-based threshold estimation,
- high-mobility fraction analysis,
- oxidative-stress-induced mobility redistribution analysis.

Generated outputs include:
- velocity link tables,
- velocity thresholds,
- per-cell velocity summaries,
- velocity-distribution figures.

---

### `04_mobility_state_analysis.ipynb`

Mobility-state classification using threshold-based and mixture-model approaches.

This notebook performs:
- slow/fast mobility classification,
- Gaussian mixture model (GMM) inference,
- mobility-state fraction analysis,
- per-cell state summaries,
- exploratory mobility-state visualization.

Generated outputs include:
- mobility-state tables,
- GMM classifications,
- state-fraction figures.

---

### `05_hmm_inference.ipynb`

Hidden Markov Model (HMM) analysis of trajectory dynamics.

This notebook performs:
- preparation of trajectory sequences,
- Gaussian HMM fitting,
- hidden-state inference,
- state-transition analysis,
- transition matrix estimation,
- per-cell HMM-state quantification.

Generated outputs include:
- HMM state assignments,
- transition matrices,
- HMM mobility-state summaries,
- state-transition figures.

---

### `06_sensitivity_analysis.ipynb`

Robustness analysis of the computational pipeline.

This notebook performs:
- testing of different edge-filter margins,
- testing of different minimum trajectory lengths,
- threshold sensitivity analysis,
- robustness evaluation of diffusion and mobility metrics.

The goal is to determine whether observed biological trends remain stable across reasonable parameter changes.

Generated outputs include:
- sensitivity-analysis tables,
- robustness figures,
- parameter-dependence summaries.

---

### `07_final_figure_generation.ipynb`

Generation of publication-ready final figures.

This notebook performs:
- loading of processed summary tables,
- standardized plot generation,
- creation of multi-panel figures,
- export of high-resolution PNG and PDF figures.

The notebook is intended for generation of reproducible figures for:
- presentations,
- posters,
- README illustrations,
- manuscripts.

Generated outputs include:
- publication-style figures,
- combined summary panels,
- export-ready graphics.

---

## Repository structure

```text
sptPALM-analysis-of-CAZ-proteins/
│
├── README.md
├── requirements.txt
├── data/
├── notebooks/
├── src/
├── figures/
├── results/
└── docs/
```

---

## Installation

Clone repository:

```bash
git clone https://github.com/YOUR_USERNAME/sptPALM-analysis-of-CAZ-proteins.git
cd sptPALM-analysis-of-CAZ-proteins
```

Create virtual environment:

### Windows

```bash
.venv\Scripts\activate
```

### Linux/macOS

```bash
source .venv/bin/activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---


## Requirements

Main Python packages:

* numpy
* pandas
* scipy
* matplotlib
* seaborn
* scikit-learn
* statsmodels
* hmmlearn
* trackpy

See requirements.txt for full list.

---

## Example workflow

### Step 1

Export trajectories from TrackMate:

* Spots CSV
* Tracks CSV

### Step 2

Place CSV files into:

```text
data/demo_data/
```

### Step 3

Run preprocessing notebook:

```text
notebooks/01_trackmate_preprocessing.ipynb
```

### Step 4

Run downstream analysis:

* diffusion analysis,
* edge filtering,
* velocity analysis,
* plotting.

---

## Results summary

The developed pipeline identified substantial heterogeneity in SNAP25 and Syntaxin mobility distributions.

Edge-associated tracking artifacts were detected and successfully removed through dedicated filtering procedures.

Velocity distributions deviated from ideal Brownian diffusion behavior and exhibited pronounced heavy tails, indicating the presence of multiple mobility states.

Oxidative stress altered mobility characteristics of synaptic proteins, with SNAP25 showing the strongest increase in highly mobile trajectory fractions after prolonged H₂O₂ exposure.

### Edge artifacts

Strong non-biological trajectory accumulation was detected near image borders.

A dedicated edge-filtering module significantly reduced this artifact.

### Diffusion heterogeneity

Protein mobility distributions strongly deviated from ideal Brownian diffusion.

Heavy-tailed velocity distributions indicate:

* heterogeneous membrane organization,
* confined diffusion,
* transient mobility-state switching.
  
### Oxidative stress effects

SNAP25 demonstrated:

* increased high-mobility fraction,
* altered velocity distributions,
* and stronger response after prolonged oxidative stress.

Syntaxin dynamics appeared more heterogeneous and condition-dependent.

---

## Current limitations

* Initial trajectory reconstruction still depends on TrackMate
* Limited number of biological replicates
* No full trajectory-state inference yet
* Spatial clustering analysis remains preliminary

---

## Future directions

Planned pipeline extensions include:

* Full Python-based trajectory linking
* Support for additional synaptic proteins
* Stronger statistical modeling
* Automated parameter optimization

---

## Data availability

This repository contains demonstration datasets, analysis notebooks, source code, and example outputs required to reproduce the computational workflow.

Original microscopy datasets are not included because of their size and institutional storage restrictions.

The metadata provided in the repository describe the structure of the analyzed experiments. Full experimental datasets may be available upon reasonable request and subject to laboratory data-sharing policies.

---

## References
1. Manley S. et al. Nature Methods (2008)
https://doi.org/10.1038/nmeth.1176
2. Betzig E. et al. Science (2006)
https://doi.org/10.1126/science.1127344
3. Hess S.T. et al. Biophysical Journal (2006)
https://doi.org/10.1529/biophysj.106.091116
4. Ershov D. et al. TrackMate 7. Nature Methods (2022)
https://doi.org/10.1038/s41592-022-01507-1
5. Jaqaman K. et al. Nature Methods (2008)
https://doi.org/10.1038/nmeth.1237
6. Kusumi A. et al. Nature Reviews Molecular Cell Biology (2012)
https://doi.org/10.1038/nrm3466

