#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 08:57:18 2018

@author: joe
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.signal import savgol_filter


dftemp = pd.read_csv('/home/joe/Code/BioSNICAR_GO/GRIS_Minerals_Mie/GRIS_dust_K.csv')
f = interpolate.interp1d(dftemp['WL'],dftemp['k'],fill_value='extrapolate')
WLnew = np.arange(0.305,4.995,0.01)
Knew = f(WLnew)

Ksmooth = savgol_filter(Knew, 11, 3) # window size 51, polynomial order 3

df = pd.DataFrame()
df['WL'] = WLnew
df['k'] = Ksmooth

plt.figure(figsize=(6,4))
plt.plot(df['WL'],df['k']),plt.xlabel('Wavelength (microns)',fontsize=16),plt.ylabel('k',fontsize=16)
plt.xticks(fontsize=16),plt.yticks(fontsize=16)
plt.savefig('/media/joe/C38C-3252/GrisMineral_k.jpg',dpi=300)




PSD = pd.read_csv('/media/joe/C38C-3252/PSD.csv',header=None)

Weights = [];

Weights.append((PSD[PSD < 0.1 ].count()) /(len(PSD)*100))
Weights.append((PSD[(PSD > 0.1) & (PSD < 0.2) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 0.2) & (PSD < 0.3) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 0.3) & (PSD < 0.4) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 0.4) & (PSD < 0.5) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 0.5) & (PSD < 0.6) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 0.6) & (PSD < 0.7) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 0.7) & (PSD < 0.8) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 0.8) & (PSD < 0.9) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 0.9) & (PSD < 1.0) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 1.0) & (PSD < 2.0) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 2.0) & (PSD < 3.0) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 3.0) & (PSD < 4.0) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 4.0) & (PSD < 5.0) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 5.0) & (PSD < 10) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 10) & (PSD < 20) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 20) & (PSD < 30) ].count())/(len(PSD)*100))
Weights.append((PSD[(PSD > 30) & (PSD < 40) ].count())/(len(PSD)*100))

Weights = np.array(Weights)

np.savetxt('PSD_Weights.csv',Weights)
