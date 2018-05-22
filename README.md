# BioSNICAR_GO
An implementation of the BioSNICAR model using geometric optics to determine ice optical properties. 

The model is a two stream radiative transfer model that predicts the albedo of an ice column with user-defined mass mixing ratios of a range of biological and non-biological particles. The optical properties of the ice grains themselves are calculated using a parameterisation of geometric optics calculations described by van Diedenhoven (2014) and the refractive indices are from Warren and brandt (2008). The radiative transfer codes are adapted from the original SNICAR model by Flanner (2007) which used Mie scattering to generate the ice optical properties. The rationale behind providing this verion is that geometic optics enables ice grains shaped as aritrarily large hexagonal plates and columns to be simulated, whereas Mie scattering applies to small spheres. The geometric ptics approach is therefore better suited to glacier ice, whereas Mie scattering is better for snow. 

# In this Repo
## Scripts
This repository contains the codes and dataset required to run the BioSNICAR_GO model. The optical properties of ice crystals are calculated using the python script 'Geometric_Optic_Ices.py'. This populates the working directory with netCDF files containing the relevant optical properties required by BioSNICAR_GO. The main BIOSNICAR_GO function is provided in the matlab script snicar8d_geom_optics.m. This is called from the driver script BioSNICAR_GO_driver.m

## Data
In this repo netCDF files are provided for ice grains between 5000 - 20000 microns in both side-length and depth, at 1000 micron resolution. Therefore, BioSNICAR_GO can be used to simulate any ice grains with dimensions in that range. The relevant files for any ice grain dimensions can be added to the working directory by running the pythin script provided.

# How to Use
Assuming the provided range of ice grain dimensions are sufficient, clone or download this repository, update the path to the working directory in line 119 of snicar8d_geometric_optics.m and then open BioSNICAR_GO_driver.m. The user-defined variables can be input at the top of the script. optional plots are commented out at the end of the script - uncoment these to view figures for subsurface light absorption etc.

# Permissions
We provide free and open access to this software on the understanding that the authors assume no responsibility for downstream usage and, while we welcome collaboration, we are not obliged to provide training or support for users. We hope that this code and data will be useful and encourage collaboration, but we provide it without warranty or guarantee, nor the implied warranty of merchantability or fitness for any particular purpose. Collaboration requests very welcome: please contact joe.cook@sheffield.ac.uk
