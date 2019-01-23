#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 15:24:05 2018

@author: joe
"""

# This code takes either a measured MAC from a csv file (300 -750 nm) or calculates
# MAC from pigment absorption coefficients and cell mass (ng/cell). It also
# optionally provides an imaginary refrcative index calculated according to a 
# mixing model adapted from Pottier et al (2005), Dauchet et al (2015) and Cook et al (2017).

# Runnign this scripts saves csv files for MAC, k etc in the working directory
# in suitable format for loading directly into the "Algae_GO.py" code that 
# returns single scattering optical proerties for the cell and saves a netcdf
# file to the Bio-SNICAR_GO lookup library

import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

plt.figure(figsize=(10,15))
plt.title("Optical Data for Glacier Algae")

# provide pigments as a mass fraction or absolute mass per cell, cell_dm_weight in ng, density in kg/m3
def bio_optical(load_MAC = True, calc_MAC = True, calc_k = True, cell_dm_weight = 0.82, chla = 0.01, chlb = 0.00066, ppro = 0.01, psyn = 0, 
                purp = 0.068, Xw = 0.5, density= 1400, nm = 1.4, savefiles = False, 
                savefilename = "name", plot_title = "title", date = "15th July 2016", plot_figs = True):
  
    data = pd.DataFrame() # set up dataframe


    if load_MAC: # choose this option to load an empirically derived MAC from file        
        MAC = pd.read_csv('/home/joe/Desktop/CW_BioSNICAR_Experiment/Empirical_MAC.csv',header=None,names=['MAC'])
        MAC = MAC[0:4695] # subsample to appropriate resolution for snicar
        MAC = MAC[0:-1:10]
        data['MAC'] = MAC['MAC'].dropna() # drop NaNs and save to dataframe
    
    if calc_MAC or calc_k: # choose this option to calculate MAC or k theoretically
        # set up empty lists
        WL = []
        Ea1 = []
        Ea2 = []
        Ea3 = []
        Ea4 = []
        Ea5 = []
        Ea1n = []
        Ea2n = []
        Ea3n = []
        Ea4n = []
        Ea5n = []
        WatRI = []
        WatRIn = []
        k_list = []
        real_list = []
    
        # Read in wavelength dependent in vivo absorption coefficients for each pigment
        # and define wavelength range
        
        WL = np.arange(300,5000,10)
        WLmeters = [i*1e-9 for i in (WL)]
    
        with open('/home/joe/Code/chlorophyll-a.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea1.append(cellf)
    
        with open('/home/joe/Code/chlorophyll-b.csv') as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea2.append(cellf)
    
        with open('/home/joe/Code/Photoprotective_carotenoids.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea3.append(cellf)
    
        with open('/home/joe/Code/Photosynthetic_carotenoids.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea4.append(cellf)
    
        with open('/home/joe/Code/phenol_MAC.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    Ea5.append(cellf)
    
        with open('/home/joe/Code/water_RI.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    WatRI.append(cellf)
    
    
       # extend spectral data down to 300nm by padding with the value at 350nm and
       # small but non-zero value above 750 nm (avoid divide by zero errors)    
        Ea1 = [Ea1[0] for _ in range(50)] + Ea1
        Ea2 = [Ea2[0] for _ in range(50)] + Ea2
        Ea3 = [Ea3[0] for _ in range(50)] + Ea3
        Ea4 = [Ea4[0] for _ in range(50)] + Ea4
      
        # extent data with zeros at nonabsoring wavelengths to 5000 nm
        for i in np.arange(751,5000,1):
            Ea1.append(0)
            Ea2.append(0)
            Ea3.append(0)
            Ea4.append(0)
            Ea5.append(0)
    
        # downsample to match SNICAR resolution) 
        for i in np.arange(1,len(Ea1),10):
            Ea1n.append(Ea1[i])
            Ea2n.append(Ea2[i])
            Ea3n.append(Ea3[i])
            Ea4n.append(Ea4[i])
            Ea5n.append(Ea5[i])
            WatRIn.append(WatRI[i])
    
        # copy water RI over zeros at non-absorbing wavTrueelengths
        Ea1n[44:-1] = WatRIn[44:-1]
        Ea2n[44:-1] = WatRIn[44:-1]
        Ea3n[44:-1] = WatRIn[44:-1]
        Ea4n[44:-1] = WatRIn[44:-1]
        Ea5n[44:-1] = WatRIn[44:-1]
    
        
        if calc_MAC:

            # convert percentage dw pigment to actual mass of pigment per cell (mg)
 
#            chla_w = chla * cell_dm_weight*1e-6
#            chlb_w = chlb * cell_dm_weight*1e-6
#            ppro_w = ppro * cell_dm_weight*1e-6
#            psyn_w = psyn * cell_dm_weight*1e-6
#            purp_w = purp * cell_dm_weight*1e-6

            # NB If data was provided in units of mg piTruegment/cell, read in
            # directly.       
            chla_w = chla
            chlb_w = chlb
            ppro_w = ppro
            psyn_w = psyn
            purp_w = purp
            
            # Multiply mass of each pigment (mg) by absorption coefficient (m2/mg)
            # at each wavelenth giving units of m2
            EW1m = [a * chla_w for a in Ea1n]
            EW2m = [a * chlb_w for a in Ea2n]
            EW3m = [a * ppro_w for a in Ea3n]
            EW4m = [a * psyn_w for a in Ea4n]
            EW5m = [a * purp_w for a in Ea5n]
            
            # Sum all pigments (m2)
            EWWm = [sum(x) for x in zip(EW1m,EW2m,EW3m,EW4m,EW5m)] 
            MAC = [i/cell_dm_weight for i in EWWm] #normalise to cell mass (ng)
            MAC = [i*1e12 for i in MAC] # convert to kg/m2
            data['MAC'] = MAC # save to dataframe


    if calc_k:
    
        # follow Pottier 2005 / Dauchet (2015) route to imaginary RI
        # multiply pigment MAC by % dw and convert from m2/mg to m2/kg
        
                # follow Pottier 2005 / Dauchet (2015) route to imaginary RI
        # multiply pigment MAC by % dw and convert from m2/mg to m2/kg
        
#        EW1 = [a * 1e6 * chla for a in Ea1n] 
#        EW2 = [a * 1e6 * chlb for a in Ea2n]
#        EW3 = [a * 1e6 * ppro for a in Ea3n]
#        EW4 = [a * 1e6 * psyn for a in Ea4n]
#        EW5 = [a * 1e6 * purp for a in Ea5n]
        
        # Sum all pigments
#        EWW = [sum(x) for x in zip(EW1,EW2,EW3,EW4,EW5)]
        
        # if data provided in mg pigment per cell, divide by cell weight in mg
        # to get mass fraction
        
        EW1 = [a  * chla/(cell_dm_weight*1e-6) for a in Ea1n] 
        EW2 = [a  * chlb/(cell_dm_weight*1e-6) for a in Ea2n]
        EW3 = [a  * ppro/(cell_dm_weight*1e-6) for a in Ea3n]
        EW4 = [a  * psyn/(cell_dm_weight*1e-6) for a in Ea4n]
        EW5 = [a  * purp/(cell_dm_weight*1e-6) for a in Ea5n]
        
        # Sum all pigments
        EWW = [sum(x) for x in zip(EW1,EW2,EW3,EW4,EW5)]
        EWW = [i*1e6 for i in EWW] # unit conversion to m2/kg
        
        # Calculate imaginary refrcative index (k)
        for i in np.arange(0,len(WL),1):
  #          k = (((1 - Xw) / Xw) * (WLmeters[i]/np.pi*4) * density * EWW[i]) #original Pottier equation
            k = (Xw * WatRIn[i]) + ((1 - Xw) * (WLmeters[i]/np.pi*4) * density * EWW[i]) # Cook (2018) updated version
            k_list.append(k)
            real_list.append(nm)
        # append real and imaginary RI to dataframe    
        data['Imag'] = k_list
        data['Real'] = real_list
        
        
        
    if savefiles: # optional save dataframe to csv files
        data.to_csv('/home/joe/Desktop/CW_BioSNICAR_Experiment/Cell_optics_dataset.csv')
        data['Imag'].to_csv('/home/joe/Desktop/CW_BioSNICAR_Experiment/{}_KK.csv'.format(savefilename),header=None,index=False)
        data['MAC'].to_csv('/home/joe/Desktop/CW_BioSNICAR_Experiment/{}_MAC.csv'.format(savefilename),header=None,index=False)
    
    if plot_figs:

        plt.subplot(2,1,1)
        plt.plot(WL[0:220],MAC[0:220],label='{}'.format(date)),plt.xlim(300,750)
        plt.xticks(fontsize=16), plt.yticks(fontsize=16)
        plt.xlabel('Wavelength',fontsize=16),plt.ylabel('MAC (kg/m^3)',fontsize=16)
        plt.legend(loc='best')
        plt.tight_layout()

        plt.subplot(2,1,2)
        plt.plot(WL[0:220],k_list[0:220],label='{}'.format(date)), plt.xlim(300,750)
        plt.xticks(fontsize=16), plt.yticks(fontsize=16)
        plt.xlabel('Wavelength',fontsize=16),plt.ylabel('K (dimensionless)',fontsize=16)
        plt.legend(loc='best')
        plt.tight_layout()

    
        return k_list, MAC, data

# NB pigment data is provided here in units of mg per cell      
k_list, MAC, data = bio_optical(
        load_MAC= True, 
        calc_MAC = False, 
        calc_k = True, 
        cell_dm_weight=1.89, # unit = ng
        chla = 3.51E-9, 
        chlb = 4.52E-9, 
        ppro = 5.725E-9, 
        psyn = 0, 
        purp = 4E-8, 
        Xw = 0.8, 
        density= 1400, 
        nm = 1.4, 
        savefiles = True, 
        savefilename = "CW_bio_1", 
        date = "15th July 2016",
        plot_figs = True)

