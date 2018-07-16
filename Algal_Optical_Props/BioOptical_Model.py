#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 15:05:56 2018

@author: joe
"""


# Code to take empirically derived values for bio-optical parameters for algal cells
# and use them to determine the wavelength-dependent imaginary part (k) of the
# complex index of refraction n' = n + ik
# This is then fed into Lorenz-Mie equations to derive absorption cross 
# section, qqabs.

# P = density of dry material (kg m^-3)
# Ea(i,wl) = in vivo absorption coefficient for each pigment per wavelength (m^2 kg^-1)
# Wi = dry mass fraction of pigment in cell (kg kg^-1 biomass)
# Xw = water fraction in each cell (dimensionless)
# Np = cell numer density (m^-3)
# V32 = mean sauter diameter (m^3)
# Cx = biomass dry weight concentration (g/L which is equal to kg/m3)
# wlmin = minimum wavelength in measured spectrum (nm)
# wlmax = maximum wavelength in measured spectrum (nm)
# dw = spectral resolution (e.g. 1nm) 


import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
from scipy.integrate import quad
import csv
import pandas as pd
from math import exp
from miepython import mie



WL = []
WLnano=[]
Ea1 = []
Ea2 = []
Ea3 = []
Ea4 = []
Ea5 = []
EW1 = []
EW2 = []
EW3 = []
EW4 = []
EW5 = []
K=[]
klist = []
Qlist = []


# Read in wavelength dependent in vivo absorption coefficients for each pigment
# from excel documents in 'C:/documents/Python Scripts'. 
# NB multiplied by 10E-10 to convert into m

with open('wavelengths.csv')as f:
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        for cell in row:
            cellf = float(cell)
            cellf_m = cellf*10E-10
            WL.append(cellf_m)

## expand WL array down to 300nm so that arrays are all same length (<350 discarded later)
 
x = [i*10E-10 for i in range(300,350)] 
x.extend(WL)
WL = x
 

# Also create new array containing wavelengths in nanometers rather than meters
# (for integration and plotting)

WLnano = [a*1e9 for a in (WL)]


with open('chlorophyll-a.csv')as f:
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        for cell in row:
            cellf = float(cell)
            Ea1.append(cellf)


with open('chlorophyll-b.csv') as f:
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        for cell in row:
            cellf = float(cell)
            Ea2.append(cellf)

with open('Photoprotective_carotenoids.csv')as f:
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        for cell in row:
            cellf = float(cell)            
            Ea3.append(cellf)
        
    
with open('Photosynthetic_carotenoids.csv')as f:
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        for cell in row:
            cellf = float(cell)
            Ea4.append(cellf)


with open('Purpurogallin.csv')as f:
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        for cell in row:
            cellf = float(cell)*1e6
            Ea5.append(cellf)

# extend spectral data down to 300nm by padding with the value at 350nm

Ea1 = [Ea1[0] for _ in range(50)] + Ea1
Ea2 = [Ea2[0] for _ in range(50)] + Ea2
Ea3 = [Ea3[0] for _ in range(50)] + Ea3
Ea4 = [Ea4[0] for _ in range(50)] + Ea4
Ea5 = [Ea5[0] for _ in range(50)] + Ea5


# if statement below means the plotting beneath it will only be executed if 
# this script is run independently, not when called by other functions.         


plt.plot(WLnano,Ea1,label='Chlorophyll a')
plt.plot(WLnano,Ea2,label='Chlorophyll b')
plt.plot(WLnano,Ea3,label='Secpndary carotenoids')
plt.plot(WLnano,Ea4,label='Primary carotenoids')
plt.plot(WLnano,Ea5,label='Purpurogallin')
plt.xlabel('Wavelengths nm')
plt.ylabel('In vivo absorption coefficient')
plt.xlim(300,750)
plt.legend(loc='best')
plt.show()

# for full model equations see Pottier et al (2005). 
# Assign mass fraction of each pigment in cell. W1 = Chll a, W2 = chll b,
# W3 = photoprotective carotenoids, W4 = photosynthetic carotenoids, W5 = phenol

W1 = 0.01
W2 = 0.00066
W3 = 0.01
W4 = 0
W5 = 0.068

# Multiply mass fraction of each pigment by absorption coefficient at each wavelenth

EW1 = [a * W1 for a in Ea1]
EW2 = [a * W2 for a in Ea2]
EW3 = [a * W3 for a in Ea3]
EW4 = [a * W4 for a in Ea4]
EW5 = [a * W5 for a in Ea5]

# Sum all pigments
EWW = [sum(x) for x in zip(EW1,EW2,EW3,EW4,EW5)]

Xw = 0.8 # water fraction in cell (can be assumed constant 0.8)
density = 1400 # density of dry material (can be assumed constant 1400 kgm-3)
nm = 1.5 # real part of RI

# Calculate refractive index, K
for i in WL:
    k = (i/(4*np.pi)) * density * ((1-Xw)/Xw) 
    klist.append(k)
    
K = [a*b for a,b in zip(EWW,klist)]
K = np.array(K)

plt.plot(WL[50:-1],K[50:-1])


KK = pd.DataFrame()
KK['WL'] = WL
KK['KK'] = K
#KK.to_csv('Real_Cell_1.csv')
