#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 15:24:05 2018

@author: joe
"""

# This code has been adapted from the original bio-optical model presented in 
# Cook et al 2017 such that it applies a volume weighted average for the RI of
# water and a summation of the other pigments when Xw is 1 or 0 instead of going
# to infinity in the latter case. This has the advantage of giving non-zero k
# and non-unity single scattering albeod at wavelengths where the pigments are
# non-absorbing. It also means the cells look like water at those non-absorbing
# wavelengths. The real part of the refractove index is changed from 1.5 to 1.4 
# based on Dauchet et al (2015).

# Function options include values for each pigment expressed as % total dry mass of the cell:
# chla = chlorophyll a, chlb = chlorophyll b, pprp = photoprotective carotenoids, psyn = photosynthetic carotenoids,
# purp = purpurogallin-like phenolic pigment

# dm_weight = weight of dry material for one cell in units of ng (mean from CW samples = 0.82ng)

# The options m2kg or m2mg are options for scaling to either unit (m^2 per kg or m^2 per mg)

import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd


def bio_optical(load_MAC = True, calc_MAC = True, calc_k = True, dm_weight = 0.82, chla = 0.01, chlb = 0.00066, ppro = 0.01, psyn = 0, 
                purp = 0.068, Xw = 0.5, density= 1400, nm = 1.4, savefiles = False, 
                savefilename = "name", plot_title = "title", plot_figs = True):

    
    data = pd.DataFrame() # set up dataframe



    if load_MAC: # choose this option to load an empirically derived MAC from file        
        MAC = pd.read_csv('/home/joe/Desktop/Empirical_MAC.csv',header=None,names=['MAC'])
        MAC = MAC[0:4695]
        MAC = MAC[0:-1:10]
        data['MAC'] = MAC['MAC'].dropna()
    

    if calc_MAC or calc_k:
        
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
        # and define wavelength range (k calculation takes WL in meters)
    
        WL = np.arange(300,5000,10)
        WLmeters = [i*1e-9 for i in (WL)]
    
        with open('/home/joe/Code/chlorophyll-a.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    Ea1.append(cellf)
    
        with open('/home/joe/Code/chlorophyll-b.csv') as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    Ea2.append(cellf)
    
        with open('/home/joe/Code/Photoprotective_carotenoids.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    Ea3.append(cellf)
    
        with open('/home/joe/Code/Photosynthetic_carotenoids.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    Ea4.append(cellf)
    
        with open('/home/joe/Code/Purpurogallin.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e6
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
        Ea5 = [Ea5[0] for _ in range(50)] + Ea5
    
        for i in np.arange(751,5000,1):
            Ea1.append(0)
            Ea2.append(0)
            Ea3.append(0)
            Ea4.append(0)
            Ea5.append(0)
    
        # downsample to match SNICAR wavelengths)
    
        for i in np.arange(1,len(Ea1),10):
            Ea1n.append(Ea1[i])
            Ea2n.append(Ea2[i])
            Ea3n.append(Ea3[i])
            Ea4n.append(Ea4[i])
            Ea5n.append(Ea5[i])
            WatRIn.append(WatRI[i])
    
        # copy water RI over zeros at non-absorbing wavelengths
        Ea1n[44:-1] = WatRIn[44:-1]
        Ea2n[44:-1] = WatRIn[44:-1]
        Ea3n[44:-1] = WatRIn[44:-1]
        Ea4n[44:-1] = WatRIn[44:-1]
        Ea5n[44:-1] = WatRIn[44:-1]
    
        # Calculate refractive indices using model from Cook et al (2017 TC) adapted so
        # that all components always sum to 1, and additional pigment is at cost of
        # water. Water RI included and mixed proportionally according to Xw
        # include RI for water
        
        if calc_MAC:

            # convert percentage dw pigment to actual mass of pigment
            chla_w = chla * dm_weight
            chlb_w = chla * dm_weight
            ppro_w = chla * dm_weight
            psyn_w = chla * dm_weight
            purp_w = chla * dm_weight
            
            # Multiply mass of each pigment by absorption coefficient at each wavelenth
            EW1m = [a * chla_w for a in Ea1n]
            EW2m = [a * chlb_w for a in Ea2n]
            EW3m = [a * ppro_w for a in Ea3n]
            EW4m = [a * psyn_w for a in Ea4n]
            EW5m = [a * purp_w for a in Ea5n]
            
            # Sum all pigments
            EWWm = [sum(x) for x in zip(EW1m,EW2m,EW3m,EW4m,EW5m)]
            MAC = [i/0.82 for i in EWWm]
            data['MAC'] = MAC


    if calc_k:
    
        # follow Pottier 2005 / Dauchet (2015) route to imaginary RI
        EW1 = [a * chla for a in Ea1n]
        EW2 = [a * chlb for a in Ea2n]
        EW3 = [a * ppro for a in Ea3n]
        EW4 = [a * psyn for a in Ea4n]
        EW5 = [a * purp for a in Ea5n]
        
        # Sum all pigments
        EWW = [sum(x) for x in zip(EW1,EW2,EW3,EW4,EW5)]

        # Calculate imaginary refrcative index (k)
        for i in np.arange(0,len(WL),1):
#            k = (((1 - Xw) / Xw) * (WLmeters[i]/np.pi*4) * density * EWW[i]) original Pottier equation
            k = (Xw * WatRIn[i]) + ((1 - Xw) * (WLmeters[i]/np.pi*4) * density * EWW[i])
            k_list.append(k)
            real_list.append(nm)
            

        data['Imag'] = k_list
        data['Real'] = real_list
            
            
    if savefiles:
        data.to_csv('/home/joe/Desktop/CW_BioSNICAR_Experiment/Cell_optics_dataset.csv')
        data['Imag'].to_csv('/home/joe/Desktop/CW_BioSNICAR_Experiment/{}_KK.csv'.format(savefilename),header=None,index=False)
        data['MAC'].to_csv('/home/joe/Desktop/CW_BioSNICAR_Experiment/{}_MAC.csv'.format(savefilename),header=None,index=False)

    
    if plot_figs:
        plt.figure(figsize=(8,8))
        plt.plot(WL,MAC)
        plt.xticks(fontsize=16), plt.yticks(fontsize=16)
        plt.xlabel('Wavelength',fontsize=16),plt.ylabel('MAC (kg/m^3)',fontsize=16)
        plt.title(plot_title,fontsize=16)
        plt.tight_layout()

        plt.figure(figsize=(8,8))
        plt.plot(WL,k_list)
        plt.xticks(fontsize=16), plt.yticks(fontsize=16)
        plt.xlabel('Wavelength',fontsize=16),plt.ylabel('K (dimensionless)',fontsize=16)
        plt.title(plot_title,fontsize=16)
        plt.tight_layout()
    
        # Plots
        plt.figure(figsize=(8,8))
        plt.plot(WL,Ea1n,label='Chlorophyll a')
        plt.plot(WL,Ea2n,label='Chlorophyll b')
        plt.plot(WL,Ea3n,label='Secondary carotenoids')
        plt.plot(WL,Ea4n,label='Primary carotenoids')
        plt.plot(WL,Ea5n,label='Purpurogallin-phenolic pigment')
        plt.xlabel('Wavelengths nm',fontsize=16)
        plt.ylabel('In vivo mass absorption coefficient (m$^2$/kg)',fontsize=16)
            
        plt.xlim(300,750), plt.xticks(fontsize=16),plt.yticks(fontsize=16)
        plt.legend(loc='best',fontsize=16)  
    
        return KK, k_list, MAC, data


for i in np.arange(0,5,1):

    if i == 0:  
        
        KK, k_list, MAC, data = bio_optical(load_MAC= True, calc_MAC = False, 
        calc_k = True, dm_weight=0.82, chla = 0.00325, chlb = 0.00049, 
        ppro = 0.00489, psyn = 0, purp = 0.0306, Xw = 0.5, density= 1400, 
        nm = 1.4, savefiles = True, savefilename = "CW_bio_1", 
        plot_title = "Imaginary part of refractive index for algae on 15$^{th}$ July 2016", 
        plot_figs = True)
    
    if i == 1:   
        
        KK, k_list, MAC, data = bio_optical(load_MAC= True, calc_MAC = False, 
        calc_k = True, dm_weight=0.82, chla = 0.00323, chlb = 0.00048, 
        ppro = 0.00385, psyn = 0, purp = 0.0357, Xw = 0.5, density= 1400, 
        nm = 1.4, savefiles = True, savefilename = "CW_bio_2", 
        plot_title = "Imaginary part of refractive index for algae on 15$^{th}$ July 2016", 
        plot_figs = True)
                
    if i == 2:        
        
        KK, k_list, MAC, data = bio_optical(load_MAC= True, calc_MAC = False, 
        calc_k = True, dm_weight=0.82, chla = 0.0019, chlb = 0.00028, 
        ppro = 0.00227, psyn = 0, purp = 0.036, Xw = 0.5, density= 1400, 
        nm = 1.4, savefiles = True, savefilename = "CW_bio_3", 
        plot_title = "Imaginary part of refractive index for algae on 15$^{th}$ July 2016", 
        plot_figs = True)
        
    if i == 3:        
        
        KK, k_list, MAC, data = bio_optical(load_MAC= True, calc_MAC = False, 
        calc_k = True, dm_weight=0.82, chla = 0.00387, chlb = 0.00061, 
        ppro = 0.00398, psyn = 0, purp = 0.0344, Xw = 0.5, density= 1400, 
        nm = 1.4, savefiles = True, savefilename = "CW_bio_4", 
        plot_title = "Imaginary part of refractive index for algae on 15$^{th}$ July 2016", 
        plot_figs = True)
        
    if i == 4:        
        KK, k_list, MAC, data = bio_optical(load_MAC= True, calc_MAC = False, 
        calc_k = True, dm_weight=0.82, chla = 0.002, chlb = 0.00041, 
        ppro = 0.00239, psyn = 0, purp = 0.0305, Xw = 0.5, density= 1400, 
        nm = 1.4, savefiles = True, savefilename = "CW_bio_5", 
        plot_title = "Imaginary part of refractive index for algae on 15$^{th}$ July 2016", 
        plot_figs = True)
        
