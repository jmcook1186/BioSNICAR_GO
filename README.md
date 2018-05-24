# BioSNICAR_GO
An implementation of the BioSNICAR model offering the user a choice between determining ice optical properties using Mie scattering (good for fine snow grains that can be assumed spherical) or geometric optics (good for wet snow and ice). This also enables different grain shapes to be used - geometric optics for hexaginal plates and columns, Mie scattering for spheres. The original SNICAR model utilised mie scattering (Flanner et al., 2007) and an update in 2017 incorporated biological impurities (Cook et al., 2017).

The model is a two stream radiative transfer model that predicts the albedo of an snow or ice with user-defined mass mixing ratios of a range of biological and non-biological particles. In GO mode, the optical properties of the ice grains themselves are calculated using a parameterisation of geometric optics calculations described by van Diedenhoven (2014) and the refractive indices are from Warren and brandt (2008). In Mie mode the ice optical properties are calculated using Mie scatering codes using the same ic erefractive indiced. The radiative transfer codes are adapted from the original SNICAR model by Flanner (2007) which used Mie scattering to generate the ice optical properties. The rationale behind providing this verion is that geometic optics enables ice grains shaped as aritrarily large hexagonal plates and columns to be simulated, whereas Mie scattering applies to small spheres. The geometric optics approach is therefore better suited to wet snow or glacier ice, whereas Mie scattering is better for snow. 

# In this Repo
## Scripts
This repository contains the codes and dataset required to run the BioSNICAR_GO model. For geometric optics the optical properties of ice crystals are calculated using the python script 'Geometric_Optic_Ices.py'. This populates the working directory with netCDF files containing the relevant optical properties required by BioSNICAR_GO. The main BIOSNICAR_GO function is provided in the matlab script snicar8d_geom_optics.m. This is called from the driver script BioSNICAR_GO_driver.m

## Data
In this repo netCDF files are provided for hexagonal ice grains between 3000 - 20000 microns in both side-length and depth, at 1000 micron resolution, or small spherical grains between 100 - 3000 microns. Therefore, BioSNICAR_GO can be used to simulate any ice grains with dimensions in that range. The relevant files for any other ice grain dimensions can be added to the working directory by running the geometric optics or Mie scattering scripts provided here.

# How to Use
Assuming the provided range of ice grain dimensions are sufficient, clone or download this repository and then open BioSNICAR_GO_driver.m. The user-defined variables can be input at the top of the script. Optional plots are commented out at the end of the script - uncoment these to view figures for subsurface light absorption etc. The folder structure in the repository should be maintained in the local wdir, otherwise updates to path definitions in several of the scripts will be required. Also see annotations in the scripts for more details on usage and background information.

# Permissions
We provide free and open access to this software on the understanding that the authors assume no responsibility for downstream usage and, while we welcome collaboration, we are not obliged to provide training or support for users. We hope that this code and data will be useful and encourage collaboration, but we provide it without warranty or guarantee, nor the implied warranty of merchantability or fitness for any particular purpose. Collaboration requests very welcome: please contact joe.cook@sheffield.ac.uk
