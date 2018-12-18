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

# The options m2kg or m2mg are options for scaling to either unit (m^2 per kg or m^2 per mg)

import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
from scipy.integrate import quad
import csv
import pandas as pd
from math import exp



def bio_optical(chla = 0.01, chlb = 0.00066, ppro = 0.01, psyn = 0, purp = 0.068, m2kg = False, m2mg = True, savefiles = False, saveplots = True):

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
    # and defne wavelength range (k calculation takes WL in meters)

    WL = np.arange(300,5000,10)
    WLmeters = [a*1e-9 for a in (WL)]


    if m2kg == True:

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

        with open('/home/joe/Code//water_RI.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    WatRI.append(cellf)


    elif m2mg ==True:

        with open('/home/joe/Code/chlorophyll-a.csv')as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)/1e6
                    Ea1.append(cellf)

        with open('/home/joe/Code/chlorophyll-b.csv') as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)/1e6
                    Ea2.append(cellf)

        with open('/home/joe/Code/Photoprotective_carotenoids.csv')as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)/1e6
                    Ea3.append(cellf)

        with open('/home/joe/Code/Photosynthetic_carotenoids.csv')as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)/1e6
                    Ea4.append(cellf)

        with open('/home/joe/Code/Purpurogallin.csv')as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    Ea5.append(cellf)

        with open('/home/joe/Code//water_RI.csv')as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)/1e6
                    WatRI.append(cellf)


                    # extend spectral data down to 300nm by padding with the value at 350nm and
    # small but non-zero value above 750 nm (avoid divide by zero errors)

    Ea1 = [Ea1[0] for _ in range(50)] + Ea1
    Ea2 = [Ea2[0] for _ in range(50)] + Ea2
    Ea3 = [Ea3[0] for _ in range(50)] + Ea3
    Ea4 = [Ea4[0] for _ in range(50)] + Ea4
    Ea5 = [Ea5[0] for _ in range(50)] + Ea5

    for i in np.arange(751,5000,1):
        Ea1.append(0.000001)
        Ea2.append(0.000001)
        Ea3.append(0.000001)
        Ea4.append(0.000001)
        Ea5.append(0.000001)

    # downsample to match SNICAR wavelengths)

    for i in np.arange(1,len(Ea1),10):
        Ea1n.append(Ea1[i])
        Ea2n.append(Ea2[i])
        Ea3n.append(Ea3[i])
        Ea4n.append(Ea4[i])
        Ea5n.append(Ea5[i])
        WatRIn.append(WatRI[i])

    # Calculate refractive indices using model from Cook et al (2017 TC) adapted so
    # that all components always sum to 1, and additional pigment is at cost of
    # water. Water RI included and mixed proportionally according to Xw
    # include RI for water

    # Multiply mass fraction of each pigment by absorption coefficient at each wavelenth
    EW1 = [a * chla for a in Ea1n]
    EW2 = [a * chlb for a in Ea2n]
    EW3 = [a * ppro for a in Ea3n]
    EW4 = [a * psyn for a in Ea4n]
    EW5 = [a * purp for a in Ea5n]

    # Sum all pigments
    EWW = [sum(x) for x in zip(EW1,EW2,EW3,EW4,EW5)]
    Xw = 0.8 # water fraction in cell (can be assumed constant 0.8)
    density = 1400 # density of dry material (can be assumed constant 1400 kgm-3)
    nm = 1.4 # real part of RI

    # Calculate imaginary refrcative index (k)k
    for i in np.arange(0,len(WL),1):
        k = (Xw * WatRIn[i]) + ((1 - Xw) * (WLmeters[i]/np.pi*4) * density * EWW[i])
        k_list.append(k)
        real_list.append(nm)

    # Plots

    plt.figure(figsize=(8,8))
    plt.plot(WL,Ea1n,label='Chlorophyll a')
    plt.plot(WL,Ea2n,label='Chlorophyll b')
    plt.plot(WL,Ea3n,label='Secondary carotenoids')
    plt.plot(WL,Ea4n,label='Primary carotenoids')
    plt.plot(WL,Ea5n,label='Purpurogallin-phenolic pigment')
    plt.xlabel('Wavelengths nm',fontsize=16)

    if m2kg ==True:
        plt.ylabel('In vivo mass absorption coefficient (m$^2$/kg)',fontsize=16)
    elif m2mg == True:
        plt.ylabel('In vivo mass absorption coefficient (m$^2$/mg)',fontsize=16)

    plt.xlim(300,750), plt.xticks(fontsize=16),plt.yticks(fontsize=16)
    plt.legend(loc='best',fontsize=16)
    if saveplots:
        plt.savefig('/home/joe/Desktop/pigment_absorption_spectra.png')
    plt.show()
    plt.close()

    plt.figure(figsize=(8,8))
    plt.plot(WL,k_list)
    plt.xticks(fontsize=16), plt.yticks(fontsize=16)
    plt.xlabel('Wavelength',fontsize=16),plt.ylabel('K',fontsize=16)
    plt.title('Imaginary part of refractive index for algal cell',fontsize=16)
    plt.tight_layout()
    if saveplots:
        plt.savefig('/home/joe/Desktop/cell_absorption_spectrum.png')
    plt.show()

    if savefiles:
        KK = pd.DataFrame()
        KK['WL'] = WL
        KK['Imag'] = k_list
        KK['Real'] = real_list
        KK.to_csv('Real_Cell_1.csv')

    return k_list



k_list = bio_optical(m2kg = False, m2mg = True, savefiles = False, saveplots = True)