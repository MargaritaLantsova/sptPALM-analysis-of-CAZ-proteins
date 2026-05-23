# sptPALM analysis of CAZ proteins

<p align="center">
  <img src="figures/readme/trajectory.gif" width="850">
</p>

Computational pipeline for single-particle tracking PALM (sptPALM) analysis of synaptic protein mobility and membrane dynamics under different cellular perturbations.

---

## Welcome!

Hello and welcome to this repository 👋

This project focuses on computational analysis of single-molecule trajectories obtained using sptPALM microscopy. The main goal is to develop a reproducible and extensible workflow for studying nanoscale dynamics of synaptic proteins under different cellular conditions.

The repository is actively evolving together with the project itself. New analysis modules, visualization approaches, statistical methods, and trajectory-processing tools will continue to be added over time.

Current pipeline features include:
- trajectory preprocessing,
- edge-artifact filtering,
- diffusion analysis,
- velocity distribution analysis,
- mobility-state characterization,
- exploratory clustering and HMM modules.

Future updates may include:
- fully automated tracking,
- advanced state inference,
- spatial nanodomain analysis,
- and support for additional synaptic proteins.

If you are interested in:
- single-molecule microscopy,
- membrane dynamics,
- computational imaging,
- or quantitative trajectory analysis,

feel free to explore the notebooks and code.
---

## Project overview

This project implements a reproducible computational workflow for analysis of single-molecule tracking data obtained using sptPALM microscopy.

The pipeline was developed for quantitative analysis of CAZ (cytomatrix active zone) and SNARE-related proteins, including:
- SNAP25b
- Syntaxin (Syx)

under:
- control conditions,
- oxidative stress (500 μM H₂O₂),
- prolonged oxidative stress exposure (30 min).

The workflow combines:
- trajectory preprocessing,
- artifact filtering,
- diffusion analysis,
- velocity statistics,
- mobility-state characterization,
- and exploratory spatial analysis.

---

## Biological background

Synaptic membrane proteins exhibit highly heterogeneous diffusion behavior due to:
- membrane compartmentalization,
- transient confinement,
- protein-protein interactions,
- clustering,
- and dynamic reorganization of synaptic nanodomains.

Single-particle tracking PALM (sptPALM) enables nanoscale analysis of these dynamics by reconstructing trajectories of individual fluorescently labeled molecules.

However, robust analysis of sptPALM data requires careful:
- trajectory filtering,
- artifact detection,
- statistical normalization,
- and reproducible computational processing.

---

## Aim

To develop a reproducible computational pipeline for:
- quantitative analysis of sptPALM trajectories,
- diffusion coefficient estimation,
- mobility-state analysis,
- artifact detection and correction,
- and comparison of protein dynamics across experimental conditions.

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

## Main results

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

This repository contains:

* analysis scripts,
* notebooks,
* demo datasets,
* example outputs,
* and generated figures.

Raw microscopy datasets are not included due to size limitations.

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

