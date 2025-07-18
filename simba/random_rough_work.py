import sys
sys.path.append("..")

import funcs_simba as funcs
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.io.fits as fits
import fitsio


path = '../simba_files/df.pkl'
df = pd.read_pickle(path)


ssfr = df['sSFR']

ssfr = ssfr[ssfr > -np.inf]
ssfr = ssfr[ssfr < np.inf]

plt.hist(ssfr)
plt.show()