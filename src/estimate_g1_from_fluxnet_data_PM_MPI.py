#!/usr/bin/env python
"""
Fit g1 from flux data using the Penman-Monteith

NB.
- this script sub-classes the fully coupled version.
- this code uses the python MPI package to speed up fits, so expect all your
cores to go to work when you run this...

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (16.05.2016)"
__email__ = "mdekauwe@gmail.com"

import matplotlib
matplotlib.use('agg') # stop windows popping up

import os
import sys
import glob
import pandas as pd
import numpy as np
from lmfit import minimize, Parameters
import matplotlib.pyplot as plt
import calendar
import datetime as dt
from scipy.stats import pearsonr
from rmse import rmse
import scipy.stats as stats
import multiprocessing as mp
import re

from fit_medlyn_gs_model import FitMedlyn
from penman_monteith import PenmanMonteith
from estimate_pressure import estimate_pressure
from estimate_g1_from_fluxnet_data_MPI import FitFluxnetData

class FitFluxnetDataPM(FitFluxnetData): #subclass, inherits from FitFluxnetData

    def worker(self, output, flist, df_site, df_lai, processor_number):

        df_out = pd.DataFrame(columns=self.out_cols)

        for i, fname in enumerate(flist):

            d = self.get_site_info(df_site, df_lai, fname)
            print d['site'], d['yr']
            df = pd.read_csv(fname, index_col='date',
                             parse_dates={'date': ["Year","DoY","Time"]},
                             date_parser=self.date_converter)

            # files contain a rouge date from the following year, fix it.
            df = self.fix_rogue_date(df, drop=True)

            if len(df[df.GPP_f > -9000.0]) == 0:
                df_out = self.write_bad_row(df_out, d)
                continue

            (months, missing_gpp) = self.get_three_most_productive_months(df)
            if missing_gpp:
                df_out = self.write_bad_row(df_out, d)
                continue

            # kPa
            df['VPD_f'] *= self.HPA_TO_KPA

            # mol m-2 s-1
            conv = self.WM2_TO_KG_M2_S * self.KG_TO_G * self.G_TO_MOL_H20
            df["ET"] = df['LE_f'] * conv

            df = self.filter_dataframe(df, d, months)

            if 'gs_est' not in df:
                df_out = self.write_bad_row(df_out, d)
            else:
                # Filter extreme gs values which have come from extremely low VPD
                extreme = df['gs_est'].mean() + (3.0 * df['gs_est'].std())
                df = df[df['gs_est'] < extreme]

                if (len(df['gs_est']) > 50 and
                    d['pft'] != "TBD" and
                    d['lat'] != "TBD" and
                    d['pft'] != "WET"):

                    # some columns have strings, force to float otherwise we
                    # will produce a bug when we attempt to fit this bad data
                    # that will be screend anyway
                    df['VPD_f'] = df['VPD_f'].astype('float64')
                    df['CO2'] = df['CO2'].astype('float64')
                    df['GPP_f'] = df['GPP_f'].astype('float64')

                    F = FitMedlyn()
                    params = F.setup_model_params()
                    (result,
                     success) = F.minimise_params(params, df, df["gs_est"])
                    if success:
                        (d, model) = self.update_fit_stats(d, result, F, df)
                        self.make_plot(F, d, df, model)

                        # Add a tropics class
                        if self.is_tropics(d):
                            d['pft'] = "TropRF"
                            df_out = self.write_row(df_out, d,
                                                    np.mean(df['GPP_f']),
                                                    np.std(df['GPP_f']),
                                                    np.mean(df['ET']),
                                                    np.std(df['ET']))
                        else:
                            df_out = self.write_row(df_out, d,
                                                    np.mean(df['GPP_f']),
                                                    np.std(df['GPP_f']),
                                                    np.mean(df['ET']),
                                                    np.std(df['ET']))

                    else:
                        df_out = self.write_bad_row(df_out, d)
                else:
                    df_out = self.write_bad_row(df_out, d)

        output.put(df_out)


    def filter_dataframe(self, df, d, months):

        # For the fully coupled this is just filler, it gets set in the PM
        # version
        d['omega'] = -999.9
        no_G = False

        # filter three most productive months
        df = df[(df.index.month == months[0]) |
                (df.index.month == months[1]) |
                (df.index.month == months[2])]

        # save months
        d['most_prod_mths'] = str(months).strip('[]')

        # NB, I'm not slicing the df here, setting this to a copy
        df_bal = df[(df['LE_fqcOK'] == 1) &
                    (df['H_fqcOK'] == 1) &
                    (df['Rn_fqcOK'] == 1) &
                    (df['G_fqcOK'] == 1)] # soil heat flux

        top = np.sum(df_bal.LE_f + df_bal.H_f)
        bot = np.sum(df_bal.Rn_f - df_bal.G_f)
        if bot > 0.0:
            d['EBR'] = top / bot
        else:
            d['EBR'] = -999.9

        # figure out summer precip before we filter
        d['summer_precip'] = np.sum(df[df["Precip_f"]>=0.0].Precip_f)

        # filter daylight hours, good LE data, GPP, CO2
        #
        # If we have no ground heat flux, just use Rn
        if len(df[(df['G_fqcOK'] == 1)]) == 0:
            df = df[(df.index.hour >= 9) &
                (df.index.hour <= 15) &
                (df['LE_fqcOK'] == 1) &
                (df['VPD_fqcOK'] == 1) &
                (df['NEE_GPP_qcOK'] == 1) &
                (df['Precip_fqcOK'] == 1) &
                (df['WS_fqcOK'] == 1) &
                (df['ET'] > 0.01 / 1000.) & # check in mmol, but units are mol
                (df['VPD_f'] > 0.05) &
                (df['GPP_f'] > 0.0) &
                (df['Rn_fqcOK'] == 1) &
                (df['Ta_fqcOK'] == 1)]
            no_G = True
        else:
            df = df[(df.index.hour >= 9) &
                (df.index.hour <= 15) &
                (df['LE_fqcOK'] == 1) &
                (df['VPD_fqcOK'] == 1) &
                (df['NEE_GPP_qcOK'] == 1) &
                (df['Precip_fqcOK'] == 1) &
                (df['WS_fqcOK'] == 1) &
                (df['VPD_f'] > 0.05) &
                (df['GPP_f'] > 0.0) &
                (df['G_fqcOK'] == 1) &
                (df['Rn_fqcOK'] == 1) &
                (df['Ta_fqcOK'] == 1)]

        # There is an issue with missing CO2 data, which would mean that bulk
        # screening the data would unnessarily remove valid data. For the odd
        # day when CO2 is missing replace with the year mean. For cases where
        # the entire year is missing e.g. Tumbarumba (all years except 2005)
        # use the global mean
        df_good_co2 = df[df['CO2'] > 0.0]
        d['global_co2'] = self.get_global_CO2_mean(d['yr'])
        d['site_co2'] = np.mean(df_good_co2.CO2)
        if len(df_good_co2) > 0:
            replacement = d['site_co2']
        else:
            replacement = d['global_co2']
            d['site_co2'] = -999.9
        df['CO2'] = np.where(df['CO2']<0, replacement, df['CO2'])

        # Replace bad site years, i.e. where site CO2 greater or less than 15%
        # the global mean
        for key, value in self.bad_co2_sites.iteritems():
            if key == d['site'] and value == d['yr']:
                df.loc[:,"CO2"] = d['global_co2']

        #df_rain = df[df.Precip_f > 0.0]
        idx = df[df.Precip_f > 0.0].index.tolist()

        bad_dates = []
        for rain_idx in idx:
            bad_dates.append(rain_idx)
            for i in xrange(48):
                new_idx = rain_idx + dt.timedelta(minutes=30)
                bad_dates.append(new_idx)
                rain_idx = new_idx

        # There will be duplicate dates most likely so remove these.
        bad_dates = np.unique(bad_dates)

        # remove rain days...
        df = df[~df.index.isin(bad_dates)]

        # Estimate gs from inverting the penman-monteith
        df = self.penman_montieth_wrapper(d, df, no_G)

        # screen for bad data, or data I've set to bad
        df = df[(df['gs_est'] > 0.0) & (np.isnan(df['gs_est']) == False)]

        return (df)

    def penman_montieth_wrapper(self, d, df, no_G):
        # convert all units...
        vpd = df['VPD_f'] * self.KPA_TO_PA # Pa
        wind = df['WS_f'] # m s-1
        rnet = df['Rn_f'] # W m-2
        if len(df.G_f):
            G = df['G_f'] # W m-2
        else:
            G = None

        if no_G:
            G = None

        tair = df['Ta_f'] # deg C
        ustar = df['ustar'] # frictional velocity, m s-1

        if d['elev'] > -500.0: # some elevations have a dash, e.g. 450-570

            PM = PenmanMonteith(use_ustar=True)
            new_press = estimate_pressure(tair, float(d['elev']))
            lambdax = PM.calc_latent_heat_of_vapourisation(tair)

            # W m-2 -> mol m-2 s-1
            #trans = df['LE_f'] / lambdax
            conv = self.WM2_TO_KG_M2_S * self.KG_TO_G * self.G_TO_MOL_H20
            trans = df['LE_f'] * conv
            df['gs_est']  = PM.invert_penman(vpd, wind, rnet, tair,
                                             new_press, trans,
                                             ustar=ustar, G=G)

            omega = PM.calc_decoupling_coefficent(wind, tair, new_press,
                                                  df["gs_est"], ustar=ustar)
            omega = omega[(omega >= 0.0) & (omega <= 1.0)]
            if len(omega) > 1:
                d['omega'] = np.nanmean(omega)
            else:
                d['omega'] = -999.9

        else:
            df['gs_est'] = -999.9

        return df

if __name__ == "__main__":

    F = FitFluxnetDataPM(fdir="data/raw_data/LaThuile_fluxnet_data/raw_data",
                         adir="data/raw_data/LaThuile_fluxnet_data/ancillary_files/csv/raw/",
                         ofdir="data/processed/",
                         co2dir="data/raw_data/global_CO2_data/",
                         site_fname="CommonAnc_LATEST.csv",
                         lai_fname = "LAI_values.csv",
                         global_co2_fname="Global_CO2_mean_NOAA.csv",
                         ofname="g1_fluxnet_PM.csv")
    F.main()
