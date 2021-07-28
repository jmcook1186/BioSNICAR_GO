# BioSNICAR_GO

## NOTE: THIS REPOSITORY IS DEPRECATED SINCE SPRING 2020. DEVELOPMENT HAS CONTINUED IN PYTHON AT THE FOLLOWING REPOSITORY.
## www.github.com/jmcook1186/BioSNICAR_GO_PY
## I STRONGLY SUGGEST MIGRATING TO THE MOST RECENT VERSION.

An implementation of the BioSNICAR model offering the user a choice between determining ice optical properties using Mie scattering (good for fine snow grains that can be assumed spherical) or geometric optics (good for wet snow and ice). This also enables different grain shapes to be used - geometric optics for hexagonal plates and columns, Mie scattering for spheres. The original SNICAR model utilised mie scattering (Flanner et al., 2007) and an update in 2017 incorporated biological impurities (Cook et al., 2017).

The model is a two stream radiative transfer model that predicts the albedo of an snow or ice with user-defined mass mixing ratios of a range of biological and non-biological particles. In GO mode, the optical properties of the ice grains themselves are calculated using a parameterisation of geometric optics calculations described by van Diedenhoven (2014) and the refractive indices are from Warren and Brandt (2008). In Mie mode the ice optical properties are calculated using Mie scatering codes using the same ic erefractive indiced. The radiative transfer codes are adapted from the original SNICAR model by Flanner (2007) which used Mie scattering to generate the ice optical properties. The rationale behind providing this verion is that geometic optics enables ice grains shaped as aritrarily large hexagonal plates and columns to be simulated, whereas Mie scattering applies to small spheres. The geometric optics approach is therefore better suited to wet snow or glacier ice, whereas Mie scattering is better for snow. 

BioSNICAR_GO also offers the option to incorporate algae as a light absorbing impurity, importantly either in the form of snow algae and mesotaenium bergrennii (small, spherical algae approximating chloromonas, chlamydomonas, chainomonas) modelled using mie theory, or Ancylonema (large (10s - 100s microns in length), long chains of cells approximated as cylinders after Lee and Pilon, 2013) modelled using geometrical optics. The mass absorption coefficient for glacier algae has been determined empirically and the model enables the user to choose between providing an empirical MAC and empirically derived pigment abundance (as absolute mass or mass fraction) or alternatively using theoretical values. This functionality, combined with the geometric optics option for the ice matrix makes BioSNICAR_GO applicable to bare glacier ice surfaces as well as wet and dry snowpacks - a significant step forwards from the original BiOSNICAR model published in Cook et al (2017). The model currently includes options for 2 glacier algae of different user-defined dimensions and one snow algae.

# Model Structure
<img src="model_structure.jpg" width=1500>


# In this Repo
## Scripts
This repository contains the codes and datasets required to run the BioSNICAR_GO model. For geometric optics the optical properties of ice crystals are calculated using the python script 'Geometric_Optic_Ice.py'. The glacier algae optical properties are calculated using the python script 'Algae_GO.py'. This populates the working directory with netCDF files containing the relevant optical properties required by BioSNICAR_GO. The main BIOSNICAR_GO function is provided in the matlab script snicar8d_geom_optics.m and snicar_8d_mie.m. These are called from the driver script BioSNICAR_GO_driver.m.

## Data
In this repo netCDF files are provided for hexagonal ice grains between 1000 - 20000 microns in both side-length and depth, at 500 - 1000 micron resolution, or small spherical grains between 100 - 3000 microns. Therefore, BioSNICAR_GO can be used to simulate any ice grains with dimensions in that range. The relevant files for any other ice grain dimensions can be added to the working directory by running the geometric optics or Mie scattering scripts provided here. We also provide an optical property library for glacier algae of a range of dimensions. Data is also provided for a packaging effect correction applied to the empirical MAC, our MAC values for individual pigments and entire cells.

# How to Use
Assuming the provided range of ice grain dimensions are sufficient, clone or download this repository and then open BioSNICAR_GO_driver.m. The user-defined variables can be input at the top of the script. Optional plots are commented out at the end of the script - uncomment these to view figures for subsurface light absorption etc. Running the driver script will call all other relevant functions and plot/print the spectral albedo to the workspace. To add new ice or algal optics to the libraries use the Bio_Optical_Model.py and Algae_GO.py scripts. The folder structure in the repository should be maintained in the local wdir, otherwise updates to path definitions in several of the scripts will be required. Also see annotations in the scripts for more details on usage and background information.

# Permissions
This model is in active development and is changing rapidly. We provide free and open access to this software on the understanding that the authors assume no responsibility for downstream usage and, while we welcome collaboration, we are not obliged to provide training or support for users. We hope that this code and data will be useful and encourage collaboration, but we provide it without warranty or guarantee, nor the implied warranty of merchantability or fitness for any particular purpose. Any usage should cite the v1.0 release (doi: 10.5281/zenodo.2598041) and the paper: 

Cook, J. M., Tedstone, A. J., Williamson, C., McCutcheon, J., Hodson, A. J., Dayal, A., Skiles, M., Hofer, S., Bryant, R., McAree, O., McGonigle, A., Ryan, J., Anesio, A. M., Irvine-Fynn, T. D. L., Hubbard, A., Hanna, E., Flanner, M., Mayanna, S., Benning, L. G., van As, D., Yallop, M., McQuaid, J., Gribbin, T., and Tranter, M.: Glacier algae accelerate melt rates on the western Greenland Ice Sheet, The Cryosphere Discuss., https://doi.org/10.5194/tc-2019-58, in review, 2019.
 
Collaboration requests very welcome: please contact joe.cook@sheffield.ac.uk

# Development
The main scripts for BioSNICAR_GO are written for Matlab because the early development began with adaptations to some legacy Matlab code. I am currently working on translating the entire code base into Python to overcome the reliance on proprietary software.
