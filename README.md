# C-3Pol_bot

**C-3Pol_bot** is a MATLAB-based toolkit for the generation of synthetic circular Dipole Spread Function (DSF) data for high-numerical-aperture (high-NA) imaging systems and for the estimation of dipole localization and polarization parameters.

The toolkit is based on the theory and methods presented in:

> I. Herrera, M. A. Alonso, and S. Brasselet, *3D Stokes Polarimetric Imaging at Nanoscales*, arXiv preprint, arXiv:2511.11222.
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

The script saves the file `I_Cf1.mat`, which contains a `15 × 30 × N_rep` array. Each `15 × 30` image corresponds to a synthetic DSF for a single polarization state. The `N_rep` dimension represents independent realizations generated with different Poisson noise distributions.

## 2. Data Analysis.


# License

This project is licensed under the BSD 3-Clause "New" or "Revised" License.


