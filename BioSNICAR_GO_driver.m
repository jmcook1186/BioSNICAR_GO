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
coszen   = 0.50;     % if DIRECT give cosine of solar zenith angle 

% THICKNESSES OF EACH VERTICAL LAYER(array) (units: meters):
dz       = [0.003 0.02 0.02 0.02 0.02];
nbr_lyr  = length(dz);  % number of snow layers

% REFLECTANCE OF SURFACE UNDERLYING SNOW:
%   Value is applied to all wavelengths.
R_sfc    = 0.15;

% DENSITY OF EACH VERTICAL LAYER (units: kg/m3)
rho_snw(1:nbr_lyr) = [500, 500, 600, 700, 700]; 

% CHOOSE METHOD FOR DETERMINING OPTICAL PROPERTIES OF ICE GRAINS
% for small spheres choose Mie, for hexagonal plates or columns of any
% size, choose GeometricOptics
Mie = 0;
GeometricOptics = 1;

%SET ICE GRAIN DIMENSIONS
% if using Mie optical properties, set rds_snw
rds_snw = [1000,2000,3000,3000,3000];

% if using GeometricOptics, set side_length and depth
side_length(1:nbr_lyr) = [3000,5000,5000,8000,10000]; 
depth(1:nbr_lyr) = [3000,5000,5000,8000,10000];

% TOTAL NUMBER OF AEROSOL SPECIES IN MODEL
nbr_aer = 32;

% CHOOSE GLACIER ALGAE DIMENSIONS
algae_r = 5; % algae radius
algae_l = 25; % algae length
wrkdir2 = '/home/joe/Code/BioSNICAR_GO/Algal_Optical_Props/'; % working directory

stb1 = 'algae_geom_'; %name stub 1
stb2 = '.nc';  % file extansion
ancyl = strcat(wrkdir2,stb1,num2str(algae_r),'_',num2str(algae_l),stb2) % create filename string


% LOOP FOR LAP MASS MIXING RATIOS IN ICE
for x = [0] 
    
% PARTICLE MASS MIXING RATIOS (units: ng(species)/g(ice), or ppb)
% add mixing ratio of each particle per vertical layer or add 'x' to 
% loop through values defined above

mss_cnc_sot1(1:nbr_lyr)  = [0,0,0,0,0];  % uncoated black carbon
mss_cnc_sot2(1:nbr_lyr)  = [0,0,0,0,0];    % coated black carbon
mss_cnc_dst1(1:nbr_lyr)  = [0,0,0,0,0];    % dust species 1
mss_cnc_dst2(1:nbr_lyr)  = [0,0,0,0,0];    % dust species 2
mss_cnc_dst3(1:nbr_lyr)  = [0,0,0,0,0];    % dust species 3
mss_cnc_dst4(1:nbr_lyr)  = [0,0,0,0,0];    % dust species 4

mss_cnc_GRISdst1(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 1 (10 micron)
mss_cnc_GRISdst2(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 2 (20 micron)
mss_cnc_GRISdst3(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 3 (30 micron)
mss_cnc_GRISdst4(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 4 (40 micron)
mss_cnc_GRISdst5(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 5 (50 micron)
mss_cnc_GRISdst6(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 6 (60 micron)
mss_cnc_GRISdst7(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 7 (70 micron)
mss_cnc_GRISdst8(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 8 (80 micron)
mss_cnc_GRISdst9(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 9 (90 micron)
mss_cnc_GRISdst10(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 10 (100 micron)
mss_cnc_GRISdst11(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 11 (110 micron)
mss_cnc_GRISdst12(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 12 (120 micron)
mss_cnc_GRISdst13(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 13 (130 micron)
mss_cnc_GRISdst14(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 14 (140 micron)
mss_cnc_GRISdst15(1:nbr_lyr)  = [0,0,0,0,0];    % GRIS dust species 15 (150 micron)

mss_cnc_ash1(1:nbr_lyr)  = [0,0,0,0,0];    % volcanic ash species 1
mss_cnc_bio1(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 1
mss_cnc_bio2(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 2
mss_cnc_bio3(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 3
mss_cnc_bio4(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 4
mss_cnc_bio5(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 5
mss_cnc_bio6(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 6
mss_cnc_bio7(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 7
mss_cnc_ancyl(1:nbr_lyr) = [1067140,0,0,0,0]; % Realistic Cell (measured pigments, 20 micron diameter)
mss_cnc_hematite(1:nbr_lyr) = [0,0,0,0,0];   % hematite
mss_cnc_mixed_sand(1:nbr_lyr) = [0,0,0,0,0];  % mixed sand

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILE NAMES CONTAINING MIE PARAMETERS FOR ALL AEROSOL SPECIES:
fl_sot1  = 'mie_sot_ChC90_dns_1317.nc';
fl_sot2  = 'miecot_slfsot_ChC90_dns_1317.nc';
fl_dst1  = 'aer_dst_bln_20060904_01.nc';
fl_dst2  = 'aer_dst_bln_20060904_02.nc';
fl_dst3  = 'aer_dst_bln_20060904_03.nc';
fl_dst4  = 'aer_dst_bln_20060904_04.nc';

fl_GRISdst1 = 'GRISdust_10.nc';
fl_GRISdst2 = 'GRISdust_20.nc';
fl_GRISdst3 = 'GRISdust_30.nc';
fl_GRISdst4 = 'GRISdust_40.nc';
fl_GRISdst5 = 'GRISdust_50.nc';
fl_GRISdst6 = 'GRISdust_60.nc';
fl_GRISdst7 = 'GRISdust_70.nc';
fl_GRISdst8 = 'GRISdust_80.nc';
fl_GRISdst9 = 'GRISdust_90.nc';
fl_GRISdst10 = 'GRISdust_100.nc';
fl_GRISdst11 = 'GRISdust_110.nc';
fl_GRISdst12 = 'GRISdust_120.nc';
fl_GRISdst13 = 'GRISdust_130.nc';
fl_GRISdst14 = 'GRISdust_140.nc';
fl_GRISdst15 = 'GRISdust_150.nc';

fl_ash1  = 'volc_ash_mtsthelens_20081011.nc';
fl_bio1  = 'biological_1.nc'; % Biological impurity 1 (30um diameter, 1.5%chll a,10% each 1 & 2 carotenoids) )
fl_bio2  = 'biological_2.nc'; % Biological impurity 2 (30um diameter, 1.5%chll a, 5% each 1 % 2 carotenoids)
fl_bio3  = 'biological_3.nc'; % Biological impurity 3 (30um diameter, 1.5%chll a, 1% each 1 % 2 carotenoids)
fl_bio4  = 'biological_4.nc'; % Biological impurity 4 (30um diameter, 1.5%chll a only)
fl_bio5  = 'biological_5.nc'; % Biological impurity 5 (10um diameter, pigs as per bio2)
fl_bio6  = 'biological_6.nc'; % Biological impurity 6 (50um diameter, pigs as per bio2)
fl_bio7  = 'biological_7.nc'; % Biological impurity 7 (20um diameter, pigs as per bio2)
fl_ancyl = ancyl; % Biological impurity with measured pigments (inc purpurogallin), 20 micron diam
fl_hematite  = 'Hematite.nc'; % Hematite particles
fl_mixed_sand  = 'Mixed_Sand.nc'; % Mixed sand (quartz and clays)


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
            mss_cnc_dst3, mss_cnc_dst4, mss_cnc_GRISdst1, mss_cnc_GRISdst2,mss_cnc_GRISdst3,mss_cnc_GRISdst4,mss_cnc_GRISdst5,mss_cnc_GRISdst6,mss_cnc_GRISdst7,mss_cnc_GRISdst8, ...
            mss_cnc_GRISdst9,mss_cnc_GRISdst10,mss_cnc_GRISdst11,mss_cnc_GRISdst12,mss_cnc_GRISdst13,mss_cnc_GRISdst14,mss_cnc_GRISdst15,...
            mss_cnc_ash1, mss_cnc_bio1, mss_cnc_bio2,mss_cnc_bio3,mss_cnc_bio4,mss_cnc_bio5, mss_cnc_bio6, mss_cnc_bio7, mss_cnc_ancyl,mss_cnc_hematite, mss_cnc_mixed_sand, fl_sot1, ...
            fl_sot2, fl_dst1, fl_dst2, fl_dst3, fl_dst4,fl_GRISdst1,fl_GRISdst2,fl_GRISdst3,fl_GRISdst4,fl_GRISdst5,fl_GRISdst6,fl_GRISdst7,fl_GRISdst8,fl_GRISdst9,fl_GRISdst10,fl_GRISdst11,fl_GRISdst12,fl_GRISdst13,...
            fl_GRISdst14,fl_GRISdst15, fl_ash1, fl_bio1,fl_bio2,fl_bio3,fl_bio4,fl_bio5,fl_bio6, fl_bio7, fl_ancyl, fl_hematite, fl_mixed_sand);
    
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
            mss_cnc_dst3, mss_cnc_dst4, mss_cnc_GRISdst1, mss_cnc_GRISdst2,mss_cnc_GRISdst3,mss_cnc_GRISdst4,mss_cnc_GRISdst5,mss_cnc_GRISdst6,mss_cnc_GRISdst7,mss_cnc_GRISdst8, ...
            mss_cnc_GRISdst9,mss_cnc_GRISdst10,mss_cnc_GRISdst11,mss_cnc_GRISdst12,mss_cnc_GRISdst13,mss_cnc_GRISdst14,mss_cnc_GRISdst15,...
            mss_cnc_ash1, mss_cnc_bio1, mss_cnc_bio2,mss_cnc_bio3,mss_cnc_bio4,mss_cnc_bio5, mss_cnc_bio6, mss_cnc_bio7, mss_cnc_ancyl,mss_cnc_hematite, mss_cnc_mixed_sand, fl_sot1, ...
            fl_sot2, fl_dst1, fl_dst2, fl_dst3, fl_dst4,fl_GRISdst1,fl_GRISdst2,fl_GRISdst3,fl_GRISdst4,fl_GRISdst5,fl_GRISdst6,fl_GRISdst7,fl_GRISdst8,fl_GRISdst9,fl_GRISdst10,fl_GRISdst11,fl_GRISdst12,fl_GRISdst13,...
            fl_GRISdst14,fl_GRISdst15, fl_ash1, fl_bio1,fl_bio2,fl_bio3,fl_bio4,fl_bio5,fl_bio6, fl_bio7, fl_ancyl, fl_hematite, fl_mixed_sand);
    end

    % process input data:   
    wvl         = data_in(:,1);   % wavelength grid
    albedo      = data_in(:,2);   % spectral albedo
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

    % figure(2);
    % plot(wvl, sub1, 'DisplayName','0 - 1 cm');
    % plot(wvl, sub2','DisplayName','1 - 3 cm');
    % plot(wvl, sub3,'DisplayName','3 - 5 cm');
    % plot(wvl, sub4,'DisplayName','5 - 7 cm');
    % plot(wvl, sub5,'DisplayName','7 - 9 cm');
    % xlabel('Wavelength (\mum)','fontsize',20);
    % ylabel('Planar Intensity Wm-2','fontsize',20);
    % xlim([0.3 2.5]);
    % ylim([0 0.03]);
    % set(gca,'xtick',0:0.1:5,'fontsize',16);
    % set(gca,'ytick',0:0.01:0.1,'fontsize',16);
    % legend;
    % hold on

    %plot total energy intensity against depth

    % figure(3);
    % plot(depths,sub_tot_line);
    % xlabel('Depth beneath surface (m)','fontsize',20);
    % ylabel('planar intensity (Wm-2)','fontsize',20);
    % set(gca,'xtick',0:0.01:0.1,'fontsize',16);
    % set(gca,'ytick',0:0.1:1,'fontsize',16);
    % hold on

    %Report albedo
    "********* REPORTING BROADBAND ALBEDDO ******* "
    
    alb_slr % albedo over solar spectrum
    flx_abs_snw; % absorbed energy in the snowpack
    heat_rt; % radiative heating rate in K/hr

    snow_depth = sum(dz); % depth of snowpack incorporating all layers
    temp_grad = ((heat_rt(1) - heat_rt(end)))/snow_depth; % temperature gradient through snowpack

end
end
