% This code takes compelx refractive index data from Zhang et al. (2015)
% for hematite (Zhang et al.'s LQ1988 hematite) and upsamples to 1nm
% spectral resolution and trims to approrpiate size for inputting into
% snicar.

% To interface with snicar, this needs to be fed into miedriver.m to obtain
% mie optical properties. These then need to be added to a NetCDF file and
% saved in the working directory. SNICAR then needs to be updated to read
% in that NetCDF file.

% Zhang et al. (2015):
% https://www.atmos-chem-phys.net/15/12159/2015/acp-15-12159-2015.pdf


% Imaginary part of RI

WL_k = [
250
300
350
400
450
500
550
600
650
700
750
800
900
1000
1250
1500
1750
2000
2500
3000
4000
5000

];

k = [0.03
0.025
0.017
0.013
0.0085
0.0078
0.0055
0.0045
0.0045
0.004
0.004
0.004
0.004
0.004
0.005
0.0057
0.0064
0.0076
0.014
0.039
0.0067
0.018
];

WL_k = WL_k(:); % ensure column vector format
k = k(:); % ensure column vector format
xi = [300:1:5000]; % set new x values within range of original WL
xi = xi(:) % ensure column vector format

ki = interp1q(WL_k,k,xi); % interpolate to new x values
ki = ki(1:4696);
ki = smooth(ki,200);
ki = ki(:);

% Real part of RI

WL_n =[

250
300
350
400
450
500
550
600
650
700
750
800
900
1000
1250
1500
1750
2000
2500
3000
3500
4000
5000


];

n = [
    1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.53
1.52
1.52
1.52
1.52
1.52
];

WL_n = WL_n(:); % ensure column vector format
n = n(:); % ensure column vector format
xii = [300:1:5000]; % set new x values within range of original WL
xii = xii(:); % ensure column vector format

ni = interp1q(WL_n,n,xii); % interpolate to new x values
ni = ni(1:4696);
ni = ni(:);

% plot real and imaginary parts of RI
figure
plot(ki)
hold on
plot(ni)