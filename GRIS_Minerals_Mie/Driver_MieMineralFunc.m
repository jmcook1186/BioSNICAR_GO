counter = 0;

for d = 0.01:0.01:10
    
    counter = counter+1;
    
    WL = 0.305:0.01:5;

    %define cell dimensions

    CellDensity = 1400; % density 1400 kg m-3
    [ExtXCmass ssa extinction scattering absorption asymmetry] = MieMineralFunc(WL, d, CellDensity,KK);
    
    ExtXCmass_out(:,counter) = ExtXCmass;
    extinction_out(:,counter) = extinction;
    scattering_out(:,counter) = scattering;
    absorption_out(:,counter) = absorption;
    ssa_out(:,counter) = ssa;
    asymmetry_out(:,counter) = asymmetry;
    
end


csvwrite('/home/joe/Desktop/ExtXCmass.csv',ExtXCmass_out)
csvwrite('/home/joe/Desktop/extinction.csv',extinction_out)
csvwrite('/home/joe/Desktop/scattering.csv',scattering_out)
csvwrite('/home/joe/Desktop/absorption.csv',absorption_out)
csvwrite('/home/joe/Desktop/ssa.csv',ssa_out)
csvwrite('/home/joe/Desktop/assymetry.csv',asymmetry_out)



% bins
%define size bins for aveaging single scattering properties

bin1 = 1:100;
bin2 = 101:200;
bin3 = 201:300;
bin4 = 301:400;
bin5 = 401:500;
bin6 = 501:600;
bin7 = 601:700;
bin8 = 701:800;
bin9 = 801:900;
bin10 = 901:1000;

for i = 0:1:10
    var = strcat('bin',num2str(i))
    ssa_bin_out(:,i) = mean(ssa_out(var))
end

