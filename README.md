# C-3Pol_bot
C-3Pol_bot is a MATLAB-based toolkit designed for the generation of synthetic circular Dipole Spread Function (DSF) data for high-numerical-aperture (high-NA) imaging systems and for the estimation of dipole localization and polarization parameters.

# Requirements
MATLAB (All codes were tested using MATLAB R2024b)

Image Processing Toolbox

Parallel Computing Toolbox 

# Workflow 

## 1. Synthetic Data Generation

### 1.1 Compute the SGM Basis

In the `SGM_bases_computation` folder, execute `SGM_bases_computation.m`. This script generates either a single SGM basis or a set of SGM bases.

The output file `MBasea.mat` is required for the generation of synthetic data.

### 1.2 Generate Circular Dipole Spread Functions (DSFs)

In the `Synthetic Data Generation` folder, execute `ellipse_DSF.m`.

This script generates synthetic circular Dipole Spread Functions (DSFs). It simulates DSFs for a specified number of signal and background photons and incorporates Poisson noise to model photon-counting statistics.

Before running this script, ensure that the file `MBasea.mat` is present in the `Synthetic Data Generation` folder.

# License

This project is licensed under the BSD 3-Clause "New" or "Revised" License.


