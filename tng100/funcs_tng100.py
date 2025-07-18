import sys
sys.path.append("..")

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import binom
import fitsio

red_hex = '#d62728'
gray_hex = '#7f7f7f'
blue_hex = '#1f77b4'


def bin_galaxies(property1,df): # Input exact column name as in Dataframe (str)
 
  '''Returns list of arrrays containing data-AGN-fraction as a function of required property, 
  values of galaxy the property, model-AGN-fraction as a function of required property
  and the errorbars on the model.

  Inputs column name of property, manga dataframe, model dataframe,'''

  if property1 == 'sSFR':
    nbins = [-np.inf,-11.5,-8]

  else:
    nbins = [9.75,10.25,10.75,11.25,11.75]

  sorted_byprop1 = df.sort_values(by = [property1])
  binned_byprop1 = pd.cut(sorted_byprop1[property1], bins= nbins)
  df[property1 + "_bin"] = binned_byprop1

  return df


def get_eddscaled_mdot(log_bhm, log_bhmdot):

  bhmdot_in_si = 10**log_bhmdot*(1.989e30/3.156e7)
  bhm_in_si = 10**log_bhm*(1.989e30)
  G = 6.67e-11
  mp = 1.673e-27
  epsilon_r = 0.2
  sigma_t = 6.65e-29
  c = 3e8

  numer = 4*np.pi*G*bhm_in_si*mp
  denom = epsilon_r*sigma_t*c

  mdot_edd_si = numer/denom

  edd_scaled_mdot = bhmdot_in_si/mdot_edd_si

  return edd_scaled_mdot


def get_edd_ratio(log_bhm, log_bhmdot, Mdot):

  epsilon_f = 0.1

  Mdot = get_eddscaled_mdot(log_bhm, log_bhmdot)
  edd_ratio = epsilon_f*Mdot

  return edd_ratio


def get_binom_ebars(frac, size):


  ppf_low = binom.ppf(0.16,size, frac)
  ppf_median = binom.ppf(0.5,size, frac)
  ppf_high = binom.ppf(0.84,size, frac)

  frac_low = ppf_low/size
  frac_high = ppf_high/size

  ebar_low = frac - frac_low
  ebar_high = frac_high - frac

  return [ebar_low, ebar_high]



def get_Fagn_of_lcrit(df,lambda_crit):

  is_strong_ragn = []
  for i in range(len(df)):

    if df['edd_ratio'][i] > lambda_crit:
      is_strong_ragn.append(True)

    else:
      is_strong_ragn.append(False)


  df['is_strong_agn'] = is_strong_ragn


  df_mass_bins = bin_galaxies("log_mass", df)
  grped_df_mass = df_mass_bins.groupby('log_mass_bin')

  mass_center = []
  size_q = []
  size_sf = []

  frac_q = []
  frac_sf = []

  ebars_q_low = []
  ebars_q_high = []

  ebars_sf_low = []
  ebars_sf_high = []


  for mass, group in grped_df_mass:

    string = str(mass)
    split = string.split(',')
    left = float(split[0][1:])
    right = float(split[1][:-1])
    center = (left + right)/2
    mass_center.append(center)

    s_q = 0
    n_agn_q = 0

    s_sf = 0
    n_agn_sf = 0

    for i in range(len(group)):
        
      if group['log_ssfr'].iloc[i] < -11.5:

        s_q = s_q + 1

        if group['is_strong_agn'].iloc[i] == True:
          n_agn_q = n_agn_q + 1


      if group['log_ssfr'].iloc[i] > -11.5:

        s_sf = s_sf + 1

        if group['is_strong_agn'].iloc[i] == True:
          n_agn_sf = n_agn_sf + 1

    size_q.append(s_q)
    frac_q.append(n_agn_q/s_q)

    [ebar_q_low, ebar_q_high] = get_binom_ebars(n_agn_q/s_q,s_q)
    ebars_q_low.append(ebar_q_low)
    ebars_q_high.append(ebar_q_high)


    size_sf.append(s_sf)
    frac_sf.append(n_agn_sf/s_sf)

    [ebar_sf_low, ebar_sf_high] = get_binom_ebars(n_agn_sf/s_sf,s_sf)
    ebars_sf_low.append(ebar_sf_low)
    ebars_sf_high.append(ebar_sf_high)

  mass_center = np.array(mass_center)

  size_q = np.array(size_q)
  frac_q = np.array(frac_q)
  ebars_q_low = np.array(ebars_q_low)
  ebars_q_high = np.array(ebars_q_high)

  size_sf = np.array(size_sf)
  frac_sf = np.array(frac_sf)
  ebars_sf_low = np.array(ebars_sf_low)
  ebars_sf_high = np.array(ebars_sf_high)

  return [mass_center,frac_q,ebars_q_low,ebars_q_high,frac_sf,ebars_sf_low,ebars_sf_high]



def plot_nlagn_constraint(ax):

  infile = '../jiyan-schechter-emcee-edr-mssfr-kormendy13-netzer19_dust_hb-fagn-0.3.2.fits'

  Fagn_lambda_cs = fitsio.read(infile, ext='FAGN_LAMBDA_CS')
  mass_bin_centers = fitsio.read(infile, ext='MASS_BIN_CENTERS')
  ssfr_bin_centers = fitsio.read(infile, ext='SSFR_BIN_CENTERS')
  mass_bin_edges = fitsio.read(infile, ext='MASS_BIN_EDGES')
  ssfr_bin_edges = fitsio.read(infile, ext='SSFR_BIN_EDGES')
  log_Fagn_mean = fitsio.read(infile, ext='LOG_FAGN_MEAN')
  log_Fagn_std = fitsio.read(infile, ext='LOG_FAGN_STD')

  ilambda = np.where(Fagn_lambda_cs == -3.)[0][0]
  log_Fagn_mean_3 = log_Fagn_mean[:, :, ilambda]
  log_Fagn_std_3 = log_Fagn_std[:, :, ilambda]

  colors = [red_hex, gray_hex, blue_hex]
  for rissfr in np.arange(len(ssfr_bin_centers), dtype=np.int32):
      issfr = len(ssfr_bin_centers) - 1 - rissfr

      for imass in np.arange(len(mass_bin_centers), dtype=np.int32):
          if imass == 0:
              if issfr == 0:
                  label = 'Intrinsic Quiescent'
              elif issfr == 2:
                  label = 'Intrinsic Star-forming'
          else:
              label = None

          log_mean = log_Fagn_mean_3[imass, issfr]
          log_std = log_Fagn_std_3[imass, issfr]

          F_mean = 10 ** log_mean
          F_lower = 10 ** (log_mean - log_std)
          F_upper = 10 ** (log_mean + log_std)

          if issfr != 1:
            ax.fill_between(
                [mass_bin_edges[imass] + 0.02, mass_bin_edges[imass + 1] - 0.02],
                [F_lower, F_lower],
                [F_upper, F_upper],
                color=colors[issfr], alpha=0.6, label=label
          )