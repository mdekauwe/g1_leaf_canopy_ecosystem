#!/usr/bin/env python
"""
Fit g1 from flux data.

NB. this code uses the python MPI package to speed up fits, so expect all your
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

class FitFluxnetData(object):
    """
    Fit Belinda's stomatal conductance model to fluxnet data using the lmfit
    package
    """

    def __init__(self, fdir, adir, co2dir, ofdir, site_fname, lai_fname,
                 global_co2_fname, ofname):

        self.flist = glob.glob(os.path.join(fdir, "*.csv"))
        self.site_fname = os.path.join(adir, site_fname)
        self.lai_fname = os.path.join(adir, lai_fname)
        self.global_co2_fname = os.path.join(co2dir, global_co2_fname)
        self.ofname = os.path.join(ofdir, ofname)
        self.out_cols = ["Scale","g1","g1_se","r2","RMSE","n","Species",\
                         "Datacontrib","Location","latitude","longitude","PFT",\
                         "Pathway","Type","Plantform", "Leafspan","Tregion"]

        # bad CO2 sites -> see python scripts/find_bad_CO2_sites.py
        self.bad_co2_sites = {'BE-Vie': '1998', 'CN-Hny': '2005',
                              'CN-Hny': '2006', 'DK-Lva': '2005',
                              'DK-Sor': '2004', 'ES-ES1': '2003',
                              'IE-Dri': '2005', 'IT-Non': '2002',
                              'PT-Esp': '2004', 'UK-Gri': '2006',
                              'UK-Her': '2006', 'US-SRM': '2004',
                              'US-SRM': '2005', 'US-SRM': '2006',
                              'US-Wkg': '2004', 'US-Wkg': '2006',
                              'ZA-Kru': '2003'}
        self.out_cols = ["site","name","country","year","latitude",\
                         "longitude","PFT", "climate_class","g1","g1_se",\
                         "n","r2","rmse","CO2", "global_CO2","summer_precip",\
                         "summer_mu_GPP_umol_m2_s", "summer_sd_GPP_umol_m2_s",\
                         "ET_mmol_m2_s","ET_sd_mmol_m2_s", "ebr","lai",\
                         "lai_min","lai_max","most_prod_mths","omega"]

        # W/m2 = 1000 (kg/m3) * 2.45 (MJ/kg) * 10^6 (J/kg) * 1 mm/day * \
        #        (1/86400) (day/s) * (1/1000) (mm/m)
        # 2.45 * 1E6 W/m2 = kg/m2/s or mm/s
        self.WM2_TO_KG_M2_S = 1.0 / ( 2.45 * 1E6 )
        self.KG_TO_G = 1000.0
        self.MOL_TO_MMOL = 1000.0
        self.G_TO_MOL_H20 = 1.0 / 18.0
        self.HPA_TO_KPA = 0.1
        self.KPA_TO_PA = 1000.0

    def main(self):

        df_site = pd.read_csv(self.site_fname)
        df_lai = pd.read_csv(self.lai_fname)

        results = self.mpi_wrapper(df_site, df_lai)

        # merge all the processor DF fits into one big dataframe
        df = pd.concat(results, ignore_index=True)

        if os.path.exists(self.ofname):
            os.remove(self.ofname)
        df.to_csv(self.ofname, index=False)

    def mpi_wrapper(self, df_site, df_lai):

        # setup multiprocessor stuff
        num_processors = mp.cpu_count()
        chunk_size = int(np.ceil(len(self.flist) / float(num_processors)))
        pool = mp.Pool(processes=num_processors)
        queue = mp.Queue() # define an output queue

        # break up the files list equally between prcoessors, of course it won't
        # quite fit eqaully so account for this
        processes = []
        for i in xrange(num_processors):

            start = chunk_size * i
            end = chunk_size * (i + 1)
            if end > len(self.flist):
                end = len(self.flist)

            # setup a list of processes that we want to run
            p = mp.Process(target=self.worker,
                           args=(queue, self.flist[start:end], df_site,
                                 df_lai, i))
            processes.append(p)

        # Run processes
        for p in processes:
            p.start()

        # OS pipes are not infinitely long, so the queuing process can get
        # blocked when using the put function - a "deadlock". The following
        # code gets around this issue

        # Get process results from the output queue
        results = []
        while True:
            running = any(p.is_alive() for p in processes)
            while not queue.empty():
                results.append(queue.get())
            if not running:
                break

        # Exit the completed processes
        for p in processes:
            p.join()

        return results

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

            # Calculate gs using pressure estimated from elevation if
            # possible, if not we will use a standard pressure assumption

            # This is a pain, but elevation could be a float or a string
            # and we need to try and catch both possibilities
            if isinstance(d['elev'], float) and d['elev'] > -500.0:
                press = estimate_pressure(df['Ta_f'], float(d['elev']))
                df['gs_est'] = df['ET'] * (press * 0.001) / df['VPD_f']

            # some elevations have dates or a dash e.g. 450-570
            elif (isinstance(d['elev'], float) == False and
                  d['elev'] != "TBD" and
                  any(c.isalpha() for c in d['elev']) == False and
                  re.search(r"-", d['elev']) == False): # some elevations
                press = estimate_pressure(df['Ta_f'], float(d['elev']))
                df['gs_est'] = df['ET'] * (press * 0.001) / df['VPD_f']
            # Calculate gs using a standard pressure
            else:
                pressure = 101.135
                df['gs_est'] = df['ET'] * pressure / df['VPD_f']

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

    def is_tropics(self, d):
        tropics = False
        if ((d['clim_grp'] == "Tropical" and d['pft'] == "DBF") or
            (d['clim_grp'] == "Tropical" and d['pft'] == "EBF") or
            (d['clim_grp'] == "Tropical" and d['pft'] == "ENF") or
            (d['lat'] >= -23.43 and d['lat'] <= 23.43 and d['pft'] == "DBF") or
            (d['lat'] >= -23.43 and d['lat'] <= 23.43 and d['pft'] == "EBF") or
            (d['lat'] >= -23.43 and d['lat'] <= 23.43 and d['pft'] == "ENF")):
            tropics = True

        return tropics

    def update_fit_stats(self, d, result, F, df):
        d['g1'] = result.params['g1'].value
        d['g1_se'] = result.params['g1'].stderr
        d['g0'] = 0.0
        model = F.gs_model(df["VPD_f"], df["GPP_f"], df["CO2"],
                           d['g0'], d['g1'])
        d['rsq'] = (pearsonr(df["gs_est"], model)[0])**2
        d['num_pts'] = len(df["gs_est"])
        d['rmse_val'] = rmse(df["gs_est"], model)

        return (d, model)

    def get_site_info(self, df_site, df_lai, fname):

        d = {}
        s = os.path.basename(fname).split(".")[0]
        d['site'] = s
        d['yr'] = os.path.basename(fname).split(".")[1]
        d['lat'] = df_site.loc[df_site.Site_ID == s,'Latitude'].values[0]
        d['lon'] = df_site.loc[df_site.Site_ID == s,'Longitude'].values[0]
        d['pft'] = df_site.loc[df_site.Site_ID == s,'IGBP_class'].values[0]
        d['clim_cl'] = df_site.loc[df_site.Site_ID == s,
                                   'Climate_class'].values[0]
        d['clim_grp'] = df_site.loc[df_site.Site_ID == s,
                                    'Climate_group'].values[0]

        # remove commas from country tag as it messes out csv output
        name = df_site.loc[df_site.Site_ID == s,'Name'].values[0]
        d['name'] = name.replace("," ,"")
        d['country'] = df_site.loc[df_site.Site_ID == s,'Country'].values[0]
        d['elev'] = df_site.loc[df_site.Site_ID == s,'Elevation'].values[0]
        d['lai'] = df_lai.loc[df_lai.sitename == s,'LAI'].values[0]
        d['lai_max'] = df_lai.loc[df_lai.sitename == s,'LAI_MAX'].values[0]
        d['lai_min'] = df_lai.loc[df_lai.sitename == s,'LAI_MIN'].values[0]

        return (d)

    def write_row(self, df_out, d, summer_gpp, summer_gpp_sd, summer_et,
                  summer_et_sd):

        row = pd.Series([d['site'], d['name'], d['country'], d['yr'], d['lat'],
                         d['lon'], d['pft'], d['clim_cl'], d['g1'], d['g1_se'],
                         d['num_pts'], d['rsq'], d['rmse_val'], d['site_co2'],
                         d['global_co2'], d['summer_precip'], summer_gpp,
                         summer_gpp_sd, summer_et*self.MOL_TO_MMOL,
                         summer_et_sd*self.MOL_TO_MMOL, d['EBR'],
                         d['lai'], d['lai_min'], d['lai_max'],
                         d['most_prod_mths'], d['omega']], index=self.out_cols)

        result = df_out.append(row, ignore_index=True)
        return result

    def write_bad_row(self, df_out, d):

        row = pd.Series([d['site'], d['name'], d['country'], d['yr'],
                         d['lat'], d['lon'], d['pft'], d['clim_cl'],
                         -999.9, -999.9, -999.9, -999.9, -999.9, -999.9,
                         -999.9, -999.9, -999.9, -999.9, -999.9, -999.9,
                         -999.9, -999.9, -999.9, -999.9, -999.9, -999.9],
                         index=self.out_cols)

        result = df_out.append(row, ignore_index=True)
        return result

    def make_plot(self, F, d, df, model):

        fig = plt.figure(figsize=(16,6))
        fig.subplots_adjust(wspace=0.05)
        fig.subplots_adjust(hspace=0.3)
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = "sans-serif"
        plt.rcParams['font.sans-serif'] = "Helvetica"
        plt.rcParams['axes.labelsize'] = 12
        plt.rcParams['font.size'] = 12
        plt.rcParams['legend.fontsize'] = 12
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12

        almost_black = '#262626'
        # change the tick colors also to the almost black
        plt.rcParams['ytick.color'] = almost_black
        plt.rcParams['xtick.color'] = almost_black

        # change the text colors also to the almost black
        plt.rcParams['text.color'] = almost_black

        # Change the default axis colors from black to a slightly lighter black,
        # and a little thinner (0.5 instead of 1)
        plt.rcParams['axes.edgecolor'] = almost_black
        plt.rcParams['axes.labelcolor'] = almost_black

        ax1 = fig.add_subplot(241, aspect='equal')
        ax2 = fig.add_subplot(243)
        ax3 = fig.add_subplot(245)
        ax4 = fig.add_subplot(246)
        ax5 = fig.add_subplot(247)

        obs = df["gs_est"]

        if np.max(obs) > np.max(model):
            maximum = np.max(obs) * 1.1
        else:
            maximum = np.max(model) * 1.1
        one2one = np.array([0, maximum])

        ax1.plot(obs, model, "ko", alpha=0.3)
        ax1.plot(one2one, one2one, 'r-')

        ax1.set_xlim(0, maximum)
        ax1.set_ylim(0, maximum)

        ax1.set_ylabel("Predicted $g_s$ (mol m$^{-2}$ s$^{-1}$)")
        ax1.set_xlabel("Inverted $g_s$ (mol m$^{-2}$ s$^{-1}$)")

        fig.suptitle("%s:%s: $g_1$ = %.2f; r$^2$ = %.2f; N = %d" % \
                     (d['site'], d['yr'], d['g1'], d['rsq'], d['num_pts']),
                     x=0.25)

        ax2.plot(df["VPD_f"], df["LE_f"], "ko", alpha=0.3)
        ax3.plot(df["VPD_f"], obs, "ko", alpha=0.3)
        ax4.plot(df["GPP_f"], obs, "ko", alpha=0.3)
        ax5.plot(df["LE_f"], obs, "ko", alpha=0.3)

        ax2.set_ylabel("LE (W m$^{-2}$)")
        ax2.set_xlabel("VPD (kPa)")

        ax3.set_xlabel("VPD (kPa)")
        ax3.set_ylabel("Observed $g_s$ (mol m$^{-2}$ s$^{-1}$)")
        ax4.set_xlabel("GPP ($\mu$mol m$^{-2}$ s$^{-1}$)")
        ax5.set_xlabel("LE (W m$^{-2}$)")

        ax1.locator_params(nbins=6)
        ax2.locator_params(nbins=6)
        ax3.locator_params(nbins=6)
        ax4.locator_params(nbins=6)
        ax5.locator_params(nbins=6)

        plt.setp(ax4.get_yticklabels(), visible=False)
        plt.setp(ax5.get_yticklabels(), visible=False)

        fig.savefig("plots/g1_fits/%s_%s_g1_fit.pdf" % (d['site'], d['yr']),
                    bbox_inches='tight', pad_inches=0.1)
        #plt.show()
        plt.close(fig)

    def get_three_most_productive_months(self, df):
        # filter three most productive months
        df_m = df.resample("M").mean()
        missing_gpp = False

        try:
            df_m = df_m.sort_values("GPP_f", ascending=False)[:3]
            months = df_m.index.month
        except KeyError:
            missing_gpp = True
            months = None

        return (months, missing_gpp)

    def filter_dataframe(self, df, d, months):

        # For the fully coupled this is just filler, it gets set in the PM
        # version
        d['omega'] = -999.9

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
        df = df[(df.index.hour >= 9) &
                (df.index.hour <= 15) &
                (df['LE_fqcOK'] == 1) &
                (df['VPD_fqcOK'] == 1) &
                (df['NEE_GPP_qcOK'] == 1) &
                (df['Precip_fqcOK'] == 1) &
                #(df['WS_fqcOK'] == 1) &
                (df['ET'] > 0.01 / 1000.) & # ET in mol, check is in mmol
                (df['VPD_f'] > 0.05) &
                (df['GPP_f'] > 0.0)]

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

        # screen for bad data, or data I've set to bad
        df = df[(df['gs_est'] > 0.0) & (np.isnan(df['gs_est']) == False)]

        return (df)

    def get_global_CO2_mean(self, yr):

        df = pd.read_csv(self.global_co2_fname, sep=",", parse_dates=False,
                        index_col=[0])
        row = df.query("year == %s" % yr)
        return row['mean'].values[0]

    def date_converter(self, *args):
        year = int(float(args[0]))
        doy = int(float(args[1]))
        # in leap years the rogue date from the following year will be 367
        # as we are correctin this below it really doesn't matter what we set
        # it to but lets make it 366 for now so that the code doesn't barf
        if doy == 367:
            doy = 366

        hour = int(args[2].split(".")[0])
        minutes = int((float(args[2]) - hour) * 60.)
        date = "%s %s %s:%s" % (year, doy, hour, minutes)

        return pd.datetime.strptime(date, '%Y %j %H:%M')

    def fix_rogue_date(self, df, drop=False):
        files_year = np.median(df.index.year)

        if drop:
            # drop 30 min slot from following year
            df = df[df.index.year == files_year]
        else:
            # files contain a rouge date from the following year, fix it.
            dates = pd.Series([pd.to_datetime(date) \
                               for date in df.index]).values
            fixed_date = "%s %s %s:%s" % (int(files_year + 1), 1, 0, 0)
            dates[-1] = pd.datetime.strptime(fixed_date, '%Y %j %H:%M')
            df = df.reindex(dates)

        return df

    def remove_from_list(self, flist, bad_files):
        return [i for i in flist if i not in bad_files]



if __name__ == "__main__":

    F = FitFluxnetData(fdir="data/raw_data/LaThuile_fluxnet_data/raw_data",
                       adir="data/raw_data/LaThuile_fluxnet_data/ancillary_files/csv/raw/",
                       ofdir="data/processed/",
                       co2dir="data/raw_data/global_CO2_data/",
                       site_fname="CommonAnc_LATEST.csv",
                       lai_fname = "LAI_values.csv",
                       global_co2_fname="Global_CO2_mean_NOAA.csv",
                       ofname="g1_fluxnet.csv")
    F.main()
