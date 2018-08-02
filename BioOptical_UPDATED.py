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


import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
from scipy.integrate import quad
import csv
import pandas as pd
from math import exp

# Set up empty lists

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
Ea6 = []
EW1 = []
EW2 = []
EW3 = []
EW4 = []
EW5 = []
WatRI = []
WatRIn = []
k_list = []
real_list = []

# Read in wavelength dependent in vivo absorption coefficients for each pigment 
# and defne wavelength range (k calculation takes WL in meters)

WL = np.arange(300,5000,10)
WLmeters = [a*1e-9 for a in (WL)]

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

with open('/home/joe/Desktop/water_RI.csv')as f:
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

# Define pigment weighting (% of total dry mass of cell)    
W1 = 0.01
W2 = 0.00066
W3 = 0.01
W4 = 0
W5 = 0.068

# Multiply mass fraction of each pigment by absorption coefficient at each wavelenth
EW1 = [a * W1 for a in Ea1n]
EW2 = [a * W2 for a in Ea2n]
EW3 = [a * W3 for a in Ea3n]
EW4 = [a * W4 for a in Ea4n]
EW5 = [a * W5 for a in Ea5n]

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
plt.plot(WL,Ea3n,label='Secpndary carotenoids')
plt.plot(WL,Ea4n,label='Primary carotenoids')
plt.plot(WL,Ea5n,label='Purpurogallin')
plt.xlabel('Wavelengths nm')
plt.ylabel('In vivo mass absorption coefficient (m2/kg)')
plt.xlim(300,750)
plt.legend(loc='best')
plt.show()

plt.figure(figsize=(8,8))
plt.plot(WL,k_list)
plt.xlabel('Wavelength'),plt.ylabel('K')


KK = pd.DataFrame()
KK['WL'] = WL
KK['Imag'] = k_list
KK['Real'] = real_list
KK.to_csv('Real_Cell_1.csv')
