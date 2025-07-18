import sys
sys.path.append("..")

import funcs_eagle as funcs
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.io.fits as fits
import fitsio


path = '../eagle_files/df.pkl'
df = pd.read_pickle(path)


plt.hist(df['log_ssfr'][(df['log_ssfr'] > -30) & (df['log_ssfr'] < np.inf)])
plt.show()