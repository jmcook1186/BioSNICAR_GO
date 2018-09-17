% Code to take values from Mie solver and use it to update NetCDF file for
% ice, water or impurity

% Inputs: filename of destination NetCDF file required. All other variables
% should be populated and available in workspace after running Mie solver
% (i.e. mieIce.m, mieWATER.m or mie.m)
% 
% Joseph Cook, Feb 2017, University of Sheffield, UK.


filename = '/home/joe/Code/BioSNICAR_GO/Algal_Optical_Props/snw_alg_40.nc';
rr(1:200) = 0.002; %radius in m

ncwrite(filename,'rds',rr) % particle radius
%ncwrite(filename,'ext_xsc',ExtXC) % extinction cross section
%ncwrite(filename,'sca_xsc',ScaXC) % scattering cross section
%ncwrite(filename,'abs_xsc',AbsXC) % absorption cross section
ncwrite(filename,'ext_cff_mss',ExtXCmass) % mass extinction cross section
%ncwrite(filename,'sca_cff_mss',ScaXCmass) % mass scattering cross section
%ncwrite(filename,'abs_cff_mss',AbsXCmass) % mass absorption cross section
%ncwrite(filename,'ext_cff_vlm',ExtXCvol) % volume extinction cross section
%ncwrite(filename,'sca_cff_vlm',ScaXCvol) % volume scattering cross section
%ncwrite(filename,'abs_cff_vlm',AbsXCvol) % volume absorption cross section
ncwrite(filename,'ss_alb',ssa) % single scattering albedo
ncwrite(filename,'asm_prm',asymmetry) % assymetry parameter
ncwrite(filename,'rds_swa',1.00001) % surface weighted radius (analytic)
ncwrite(filename,'rds_swr',1.00001) % surface weighted radius (resolved)
ncwrite(filename,'rds_nma',1) % analytic number-mean radius
ncwrite(filename,'gsd',5e-6) % geometric SD of lognormal distribution
ncwrite(filename,'prt_dns',1) % particle density
