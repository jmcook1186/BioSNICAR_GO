% Driver routine for BioSNICAR_GO.
% Written by Joseph Cook, 2018 (joe.cook@sheffield.ac.uk)

% This code adds significant functionality to the original SNICAR model
% (Flanner et al, 2007) by enabling users to model large, aspherical ice
% grains as well as small, spherical ones as per SNICAR. The model also
% offers incorporation of biological impurities as per BioSNICAR (Cook et
% al., 2017).

% The user chooses to access ice grains whose optical properties are
% calculated using Mie scattering calculations (assuming grains are small
% and spherical - good for snow) or hexagonal plates/columns of arbitrarily
% large dimensions using geometric optics (good for wet snow/ice).

% If the user chooses to use Mie scattering, the ice grain effective radii
% should then be defined under the variable name rds_snw for each vertical
% layer. 

% if the user chooses to use geometric optics, the ice grains dimensions 
% are described using the length of one side of the hexagonal face
% (side_length) and the column depth (depth). If the aspect ratio is less
% than 1 this is a hexagnal plate, if the aspect ratio is greater than one
% it is a hexagonal column. From these user-defined dimensions, the volume,
% surface area and apotherm of the hexagonal ice grains are calculated and
% reported to the command window.

% From July 2018 BioSNICAR_GO also uses the geometrical optics
% approximation to determine the optical properties of glacier algae,
% assuming them to be cylinders of dimensions defined in this driver. The
% user simply defines the radius and length of the desired glacier algal
% cell to retrieve the optical properties from a library populated using a
% python model (Algae_GO.py). The refractive indexes are determined using
% the BiOSNICAR mixing model (updated July 2018) from values measured in
% the field on the Black and Bloom project (SW GrIS).

%%%%%%%%%%  Input parameters: %%%%%%%%%%%

% BND_TYP:      Spectral grid (=1 for 470 bands. This is the
%               only functional option in this distribution)
% DIRECT:       Direct or diffuse incident radiation (1=direct, 0=diffuse)
% APRX_TYP:     Two-Stream Approximation Type:
%                  1 = Eddington
%                  2 = Quadrature
%                  3 = Hemispheric Mean
%                  NOTE: Delta-Eddington Approximation is probably
%                  best for the visible spectrum, but can provide
%                  negative albedo in the near-IR under diffuse
%                  light conditions. Hence, the Hemispheric Mean
%                  approximation is recommended for general use.
% DELTA:        1=Use Delta approximation (Joseph, 1976), 0=don't use
% coszen:       cosine of solar zenith angle (only applies when DIRECT=1)
% R_sfc:        broadband albedo of underlying surface 
%                  (user can also define a spectrally-dependent albedo below)
% dz:           array of snow layer thicknesses [meters]. Top layer has index 1. 
%                  The length of this array defines number of snow layers
% rho_snw:      array of snow layer densities [kg m-3]. Must have same length as dz
% rds_snw:      effective radii for ice grains. Only used for Mie option
% side_length:  length of one side of hexagonal plane of ice crystal [microns]. Only used for GO option
% depth:        vertical dimension of hexagonal columnar ice crystal [microns]. Only used for GO option  
% nbr_aer:      number of aerosol species in snowpack
% mss_cnc_X:    mass mixing ratio of aerosol X
% fl_X:         name of file containing optical properties for aerosol X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%%%%%%%%%%%%%%% USER DEFNED INPUTS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RADIATIVE TRANSFER CONFIGURATION:
BND_TYP  = 1;        % 1= 470 spectral bands
DIRECT   = 0;        % 1= Direct-beam incident flux, 0= Diffuse incident flux
APRX_TYP = 1;        % 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
DELTA    = 1;        % 1= Apply Delta approximation, 0= No delta
coszen   = 0.57;     % if DIRECT give cosine of solar zenith angle 

% THICKNESSES OF EACH VERTICAL LAYER(array) (units: meters):
dz       = [0.003 0.02 0.02 0.02 0.02];
nbr_lyr  = length(dz);  % number of snow layers

% REFLECTANCE OF SURFACE UNDERLYING SNOW:
%   Value is applied to all wavelengths.
R_sfc    = 0.15;

% DENSITY OF EACH VERTICAL LAYER (units: kg/m3)
rho_snw(1:nbr_lyr) = [350, 400, 500, 700, 700]; 

% CHOOSE METHOD FOR DETERMINING OPTICAL PROPERTIES OF ICE GRAINS
% for small spheres choose Mie, for hexagonal plates or columns of any
% size, choose GeometricOptics
Mie = 0;
GeometricOptics = 1;

%SET ICE GRAIN DIMENSIONS
% if using Mie optical properties, set rds_snw
rds_snw = [400,400,400,400,400];

% if using GeometricOptics, set side_length and depth
side_length(1:nbr_lyr) = [3000,3000,4000,8000,10000]; 
depth(1:nbr_lyr) = [3000,3000,4000,8000,10000];

% TOTAL NUMBER OF AEROSOL SPECIES IN MODEL
nbr_aer = 10;

% CHOOSE GLACIER ALGAE DIMENSIONS
algae_r = 3; % algae radius
algae_l = 10; % algae length
wrkdir2 = '/home/joe/Code/BioSNICAR_GO/Algal_Optical_Props/'; % working directory

stb1 = 'algae_geom_'; %name stub 1
stb2 = '.nc';  % file extansion
ancyl = strcat(wrkdir2,stb1,num2str(algae_r),'_',num2str(algae_l),stb2) % create filename string


% CHOOSE SNOW ALGAE DIAMETER
snw_algae_r = 1; % algae diameter
wrkdir2 = '/home/joe/Code/BioSNICAR_GO/Algal_Optical_Props/'; % working directory

snw_stb1 = 'snw_alg_'; %name stub 1
snw_stb2 = '.nc';  % file extansion
snw_alg = strcat(wrkdir2,snw_stb1,num2str(snw_algae_r),snw_stb2) % create filename string


% LOOP FOR LAP MASS MIXING RATIOS IN ICE
for x = [0] 
    
% PARTICLE MASS MIXING RATIOS (units: ng(species)/g(ice), or ppb)
% add mixing ratio of each particle per vertical layer or add 'x' to 
% loop through values defined above

mss_cnc_sot1(1:nbr_lyr)  =    [0,0,0,0,0];    % uncoated black carbon
mss_cnc_sot2(1:nbr_lyr)  =    [0,0,0,0,0];    % coated black carbon
mss_cnc_dst1(1:nbr_lyr)  =    [0,0,0,0,0];    % dust species 1
mss_cnc_dst2(1:nbr_lyr)  =    [0,0,0,0,0];    % dust species 2
mss_cnc_dst3(1:nbr_lyr)  =    [0,0,0,0,0];    % dust species 3
mss_cnc_dst4(1:nbr_lyr)  =    [0,0,0,0,0];    % dust species 4
mss_cnc_ash1(1:nbr_lyr)  =    [0,0,0,0,0];    % volcanic ash species 1
mss_cnc_GRISdust(1:nbr_lyr) = [0,0,0,0,0];    % GRIS dust
mss_cnc_snw_alg(1:nbr_lyr)  = [0,0,0,0,0];    % Snow Algae (spherical, C nivalis)
mss_cnc_ancyl(1:nbr_lyr) =    [0,0,0,0,0];    % Ancylonema Nordenskioldii

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILE NAMES CONTAINING MIE PARAMETERS FOR ALL AEROSOL SPECIES:
fl_sot1  = 'mie_sot_ChC90_dns_1317.nc';
fl_sot2  = 'miecot_slfsot_ChC90_dns_1317.nc';
fl_dst1  = 'aer_dst_bln_20060904_01.nc';
fl_dst2  = 'aer_dst_bln_20060904_02.nc';
fl_dst3  = 'aer_dst_bln_20060904_03.nc';
fl_dst4  = 'aer_dst_bln_20060904_04.nc';
fl_ash1  = 'volc_ash_mtsthelens_20081011.nc';
fl_GRISdust = 'GRISdust_PSD.nc';
fl_snw_alg  = snw_alg; % snow algae (c nivalis)
fl_ancyl = ancyl; % Glacier algae (ancylonema nordenskioldii)

% Check that one method for ice optical properties is selected, if not
% raise exception

if GeometricOptics == 1 && Mie ==1 

    error = "You have selected both methods for ice grain optics - please pick one!" 

elseif GeometricOptics == 0 && Mie ==0
    error = "You have not selected a method for ice grain optics, please choose either Mie or Geometric Optics"

else
    
    % run snicar in either GO or Mie mode:
    
    if GeometricOptics == 1
        
        data_in = snicar8d_GO(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, ...
            dz, rho_snw, side_length, depth, nbr_aer, mss_cnc_sot1, ...
            mss_cnc_sot2, mss_cnc_dst1, mss_cnc_dst2, ...
            mss_cnc_dst3, mss_cnc_dst4, ...
            mss_cnc_ash1, mss_cnc_GRISdust, mss_cnc_snw_alg, mss_cnc_ancyl, fl_sot1, ...
            fl_sot2, fl_dst1, fl_dst2, fl_dst3, fl_dst4, fl_ash1, fl_GRISdust, fl_snw_alg, fl_ancyl);
    
        for i = 1:1:length(dz)
            "******** REPORTING ICE GRAIN DIMENSIONS ********"
            Volume(i) = 1.5*sqrt(3)*(side_length(i)^2)*depth(i) % ice grain volume
            Area_total(i) = 3 * side_length(i) * (sqrt(3)*side_length(i)+depth(i)*2) %total surface area 
            Area(i) = Area_total(i)/4   % projected area
            apothem(i) = (2*Area(i)) / (depth(i)*6) % apothem is distance from centre point to midpoint of a side for hexagon
            diameter(i) = 2*apothem(i) % midpoint of one side to midpoint of opposite side
        end
        
    end

    if Mie == 1
        
            data_in = snicar8d_Mie(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, ...
            dz, rho_snw, rds_snw, nbr_aer, mss_cnc_sot1, ...
            mss_cnc_sot2, mss_cnc_dst1, mss_cnc_dst2, ...
            mss_cnc_dst3, mss_cnc_dst4, ...
            mss_cnc_ash1, mss_cnc_GRISdust, mss_cnc_snw_alg, mss_cnc_ancyl, fl_sot1, ...
            fl_sot2, fl_dst1, fl_dst2, fl_dst3, fl_dst4,fl_ash1, fl_GRISdust, fl_snw_alg, fl_ancyl);
    end

    % process input data:   
    wvl         = data_in(:,1);   % wavelength grid
    albedo      = data_in(:,2);   % spectral albedo
    
    
    % do not allow albedo to drop below 0

    for i = 1:1:length(albedo)
        if albedo(i) <= 0;
            albedo(i) = 0.0000001;
        end
    end
    
    alb_slr     = data_in(1,3);   % broadband albedo (0.3-5.0um)
    alb_vis     = data_in(2,3);   % visible albedo (0.3-0.7um)
    alb_nir     = data_in(3,3);   % near-IR albedo (0.7-5.0um)
    flx_abs_snw = data_in(4,3);   % total radiative absorption by all snow layers (not including underlying substrate)
    flx_abs(1)     = data_in(6,4); % top layer solar absorption
    flx_vis_abs(1) = data_in(7,4); % top layer VIS absorption
    flx_nir_abs(1) = data_in(8,4); % top layer NIR absorption
    heat_rt = data_in([1:5],6); % heating rate per layer K/hr
    F_abs = [data_in(6,4),data_in(9,4),data_in(12,4),data_in(15,4),data_in(18,4)]; % absorbed energy per layer (W/m2)
    intensity2 = data_in([1:470],[7:11]);
    
    %albedo = smooth(albedo,0.005); % add a simple smoothing function with short period
    %separate subsurface light field into layers
    sub1 = smooth(intensity2(:,1),4);
    sub2 = smooth(intensity2(:,2),4);
    sub3 = smooth(intensity2(:,3),4);
    sub4 = smooth(intensity2(:,4),4);
    sub5 = smooth(intensity2(:,5),4);

    sub1_tot = sum(sub1);
    sub2_tot = sum(sub2);
    sub3_tot = sum(sub3);
    sub4_tot = sum(sub4);
    sub5_tot = sum(sub5);
    sub_tot_line = [sub1_tot,sub2_tot,sub3_tot,sub4_tot,sub5_tot];
    depths = [dz(1),dz(1)+dz(2),dz(1)+dz(2)+dz(3),dz(1)+dz(2)+dz(3)+dz(4), dz(1)+dz(2)+dz(3)+dz(4)+dz(5)];

    figure(1);

    % make a plot of spectrally-resolved albedo:
    plot(wvl,albedo,'linewidth',2);
    xlabel('Wavelength (\mum)','fontsize',20);
    ylabel('Albedo','fontsize',20);
    set(gca,'xtick',0:0.1:5,'fontsize',16);
    set(gca,'ytick',0:0.1:1.0,'fontsize',16);
    xlim([0.3 2.5])
    ylim([0,1])
    grid on;
    hold on;

    % plot subsurface light field

    figure(2);
    hold on
    plot(wvl, sub1, 'DisplayName','0 - 0.3 cm');
    plot(wvl, sub2','DisplayName','0.3 - 1.3 cm');
    plot(wvl, sub3,'DisplayName','1.3 - 3.3 cm');
    plot(wvl, sub4,'DisplayName','3.3 - 5.3 cm');
    plot(wvl, sub5,'DisplayName','5.3 - 7.3 cm');
    xlabel('Wavelength (\mum)','fontsize',20);
    ylabel('Planar Intensity Wm-2','fontsize',20);
    xlim([0.3 2.5]);
    ylim([0 0.03]);
    set(gca,'xtick',0:0.1:5,'fontsize',16);
    set(gca,'ytick',0:0.01:0.1,'fontsize',16);
    legend;

    %plot total energy intensity against depth

    % figure(3);
    % plot(depths,sub_tot_line);
    % xlabel('Depth beneath surface (m)','fontsize',20);
    % ylabel('planar intensity (Wm-2)','fontsize',20);
    % set(gca,'xtick',0:0.01:0.1,'fontsize',16);
    % set(gca,'ytick',0:0.1:1,'fontsize',16);
    % hold on

    %Report albedo
   
    "********* REPORTING BROADBAND ALBEDO ******* "
    
    alb_slr % albedo over solar spectrum
    flx_abs_snw; % absorbed energy in the snowpack
    heat_rt; % radiative heating rate in K/hr

    snow_depth = sum(dz); % depth of snowpack incorporating all layers
    temp_grad = ((heat_rt(1) - heat_rt(end)))/snow_depth; % temperature gradient through snowpack

end
end
