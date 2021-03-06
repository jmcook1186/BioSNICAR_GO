counter = 0;

for d = 0.01:0.01:35
    
    counter = counter+1;
    
    WL = 0.305:0.01:5;

    %define cell dimensions
    CellDensity = 4287; % density 1400 kg m-3
    [ExtXCmass ssa extinction scattering absorption asymmetry] = MieMineralFunc(WL, d, CellDensity,KK);
    
    % Each single scattering optical property has its own output table - here
    % the values for each size are added to the table as a new column

    ExtXCmass_out(:,counter) = ExtXCmass;
    extinction_out(:,counter) = extinction;
    scattering_out(:,counter) = scattering;
    absorption_out(:,counter) = absorption;
    ssa_out(:,counter) = ssa;
    asymmetry_out(:,counter) = asymmetry;
    
end


%write the optical property tables to individual csv files
csvwrite('/home/joe/Code/BioSNICAR_GO/GRIS_Minerals_Mie/ExtXCmass.csv',ExtXCmass_out)
csvwrite('/home/joe/Code/BioSNICAR_GO/GRIS_Minerals_Mie/extinction.csv',extinction_out)
csvwrite('/home/joe/Code/BioSNICAR_GO/GRIS_Minerals_Mie/scattering.csv',scattering_out)
csvwrite('/home/joe/Code/BioSNICAR_GO/GRIS_Minerals_Mie/absorption.csv',absorption_out)
csvwrite('/home/joe/Code/BioSNICAR_GO/GRIS_Minerals_Mie/ssa.csv',ssa_out)
csvwrite('/home/joe/Code/BioSNICAR_GO/GRIS_Minerals_Mie/assymetry.csv',asymmetry_out)


% bins
%define size bins for aveaging single scattering properties

bin1 = 1:10; % 0 - 0.1 micron
bin2 = 101:20; % 0.1 - 0.2 micron
bin3 = 201:30; %0.2 - 0.3 micron
bin4 = 301:40; % 0.3 - 0.4 micron
bin5 = 401:50; % 0.4 - 0.5 micron
bin6 = 501:60; % 0.5 - 0.6 micron
bin7 = 601:70; % 0.6 - 0.7 micron
bin8 = 701:80; % 0.7 - 0.8 micron
bin9 = 801:90; % 0.8 - 0.9 micron
bin10 = 901:100; % 0.9 - 1 micron
bin11 = 101:200 % 1 - 2 micron
bin12 = 201:300;% 2 - 3 micron
bin13 = 301:400; % 3 - 4 micron
bin14 = 401:500; % 4 - 5 micron
bin15 = 501:1000; % 5-10 micron
bin16 = 1001:2000; %10 - 20 micron
bin17 = 2001:3000; % 20-30 micron
bin18 = 3001:4000; % 30-40 micron


% calculate the average single scattering optical properties for each bin
for i = 1:1:18
    bin = strcat('bin',num2str(i));
    ssa_bin_out(:,i) = mean(ssa_out([1:470],[bin]),2);
    ExtXCmass_bin_out(:,i) = mean(ExtXCmass_out([1:470],[bin]),2);
    extinction_bin_out(:,i) = mean(extinction_out([1:470],[bin]),2);
    scattering_bin_out(:,i) = mean(scattering_out([1:470],[bin]),2);
    absorption_bin_out(:,i) = mean(absorption_out([1:470],[bin]),2);
    asymmetry_bin_out(:,i) = mean(asymmetry_out([1:470],[bin]),2);
end

% find the bulk refractive index by taking mean of all bins weighted by
% frequency of occurrence in PSD

Weights = csvread('/home/joe/Code/BioSNICAR_GO/GRIS_Minerals_Mie/PSD_Weights.csv');
Weights = Weights/100;

for i = 1:1:18
    ssa_w_temp(:,i) = ssa_bin_out(:,i).*Weights(i);
    ExtXCmass_w_temp(:,i) = ExtXCmass_bin_out(:,i).*Weights(i);
    extinction_w_temp(:,i) = extinction_bin_out(:,i).*Weights(i);
    absorption_w_temp(:,i) = absorption_bin_out(:,i).*Weights(i);
    asymmetry_w_temp(:,i) = asymmetry_bin_out(:,i).*Weights(i);
    scattering_w_temp(:,i) = scattering_bin_out(:,i).*Weights(i);
end

weighted_ssa = sum(ssa_w_temp,2);
weighted_ExtXCmass = sum(ExtXCmass_w_temp,2);
weighted_asymmetry = sum(asymmetry_w_temp,2);
weighted_extinction = sum(extinction_w_temp,2);
weighted_scattering = sum(scattering_w_temp,2);
weighted_absorption = sum(absorption_w_temp,2);



% check output in plots

figure

subplot(3,2,1)
plot(WL,weighted_ExtXCmass)
xlabel('Wavelength (micron)')
ylabel('Cext')

subplot(3,2,2)
plot(WL,weighted_ssa)
xlabel ('wavelength (micron)')
ylabel('ssa')

subplot(3,2,3)
plot(WL,weighted_extinction)
xlabel ('wavelength (micron)')
ylabel('extinction efficiency')

subplot(3,2,4)
plot(WL,weighted_scattering)
xlabel ('wavelength (micron)')
ylabel('scattering efficiency')

subplot(3,2,5)
plot(WL,weighted_absorption)
xlabel ('wavelength (micron)')
ylabel('absorption efficiency')

subplot(3,2,6)
plot(WL,weighted_asymmetry)
xlabel ('wavelength (micron)')
ylabel('asymmetry parameter (g)')