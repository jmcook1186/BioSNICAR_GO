% Driver software for determining the mie optical parameters for biological
% cells to feed into BioSNICAR. Calls mie.m, which should be available in
% workspace. Written by Joseph Cook (University of Sheffield, UK), Feb 2017.

%User defined inputs = CELL DIAMETER (line 12) AND ABSORPTION COEFFICIENT
%(PROVIDE BY IMPORTING RELEVANT FILE TO WORKSPACE AND NAMING 'KK') 

% OUTPUTS: extinction efficiency (extinction), backscattering efficiency
% (backscattering), absorption efficiency (absorption), asymmetry parameter
% (asymmetry), backscattering:scattering ratio (q_ratio), single scattering
% albedo (ssa),extinction cross section (ExtXC), scttering cross section
% (ScaXC), absorption cross section (AbsXC), volume extinction coefficient
% (ExtXCvol), volume scattering coefficient (ScaXCvol, volume absorption
% coefficient (AbsXCvol), mass extinction coefficient (ExtXCmass), mass
% scattering coefficient (ScaXCmass), mass absorption coefficient
% (AbsXCmass)

function [ExtXCmass ssa extinction scattering absorption asymmetry] = MieMineralFunc(WL, d, CellDensity,KK)

    %define cell dimensions
    CellVol = 4/3* pi * ((d*1e-6)/2)^3; % cell vol m3
    XSArea = pi*((d*1e-6)/2)^2; % cell XS area m2

    %Calculate Mie size parameter
    XX = pi * d ./ WL;

    %call mie.m with XX and KK values
    for i = 1:1:470
        m = complex(1.45,KK(i)); % adjust real part of the RI if necessary (default for cells is 1.5 after Dauchet et al. 2005)
        x = XX(i); %read in size parameter
        [qext qsca qabs qb asy qratio] = Mie(m,x); % send size param and imaginary RI to mie.m for each wavelength
        extinction(i)=qext;
        scattering(i) = qsca;
        absorption(i) = qabs;
        backscattering(i) = qb;
        asymmetry(i) = asy;
        q_ratio(i) = qratio;
        ssa(i) = qsca/qext;
    end

    ExtXC = (extinction);
    ScaXC = (scattering);
    AbsXC = (absorption);

    ExtXCvol = (extinction/XSArea).*CellVol;
    ScaXCvol = (scattering/XSArea).*CellVol;
    AbsXCvol = (absorption/XSArea).*CellVol;

    ExtXCmass = (((extinction*XSArea)./(CellVol*CellDensity))); % MAC in kg m-3
    ScaXCmass = (((scattering*XSArea)./(CellVol*CellDensity))); % MAC in kg m-3
    AbsXCmass = (((absorption*XSArea)./(CellVol*CellDensity))); % MAC in kg m-3
    

end
