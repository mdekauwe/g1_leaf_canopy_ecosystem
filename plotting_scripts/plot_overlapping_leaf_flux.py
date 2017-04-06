#!/usr/bin/env python

"""
Where leaf gas exchange dates match with flux data, plot gs at the two scales

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (20.04.2016)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import re
import brewer2mpl

from scipy import stats
sys.path.append('src')
from estimate_g1_from_fluxnet_data_PM_MPI import FitFluxnetData
from estimate_pressure import estimate_pressure
from fit_medlyn_gs_model import FitMedlyn

import matplotlib.pyplot as plt


def main():
    """
    Two steps:
    (i) collect all the leaf and matching flux data for each of the sites.
    (ii) make big summary plot
    """

    # store all the sites...
    df_leaf_all = pd.DataFrame()
    df_flux_all = pd.DataFrame()
    sites = []

    leaf_g1 = {}
    flux_g1 = {}

    FM = FitMedlyn(fluxnet=True)
    F = GetFluxData(fdir="data/raw_data/LaThuile_fluxnet_data/raw_data",
                    adir=("data/raw_data/LaThuile_fluxnet_data/"
                          "ancillary_files/csv/raw/"),
                    ofdir="data/processed/",
                    co2dir="data/raw_data/global_CO2_data/",
                    site_fname="CommonAnc_LATEST.csv",
                    lai_fname = "LAI_values.csv",
                    global_co2_fname="Global_CO2_mean_NOAA.csv",
                    ofname="g1_fluxnet.csv")


    leaf_dir = "data/raw_data/leafgasexchange"
    fname = "WUEdatabase_13_04_2016_fixed.csv"
    df_wue = pd.read_csv(os.path.join(leaf_dir, fname))
    df_wue = add_pft_info(df_wue)

    fdir="data/raw_data/LaThuile_fluxnet_data/raw_data"

    #############
    # Hytiala
    #############
    site = "FI-Hyy"
    year = 2006
    sites.append(site)
    locn = "Hyytiala Finland"
    df_leaf, g1 = get_leaf_data(df_wue, locn)
    df_leaf["site"] = site
    leaf_g1[site] = g1

    fname = os.path.join(fdir, "%s.%s.synth.hourly.allvars.csv" % (site, year))
    df_flux = pd.DataFrame()
    df, g1 = F.main(fname)
    df_flux = pd.concat([df_flux, df])
    flux_g1[site] = g1

    df_flux["gs_pred"] = FM.gs_model(df_flux["VPD_f"], df_flux["GPP_f"],
                                     df_flux["CO2"], 0.0, g1)
    df_flux["site"] = site

    df_leaf_all = pd.concat([df_leaf_all, df_leaf])
    df_flux_all = pd.concat([df_flux_all, df_flux])



    #############
    # Le Bray
    #############
    site = "FR-LBr"
    year = 1997
    sites.append(site)
    locn = "Le Bray_Bordeaux_France"
    df_leaf, g1 = get_leaf_data(df_wue, locn)
    df_leaf["site"] = site
    leaf_g1[site] = g1

    fname = os.path.join(fdir, "%s.%s.synth.hourly.allvars.csv" % (site, year))
    df_flux = pd.DataFrame()
    df, g1 = F.main(fname)
    df_flux = pd.concat([df_flux, df])
    flux_g1[site] = g1

    df_flux["gs_pred"] = FM.gs_model(df_flux["VPD_f"], df_flux["GPP_f"],
                                     df_flux["CO2"], 0.0, g1)
    df_flux["site"] = site

    df_leaf_all = pd.concat([df_leaf_all, df_leaf])
    df_flux_all = pd.concat([df_flux_all, df_flux])



    #############
    # Griffin
    #############
    site = "UK-Gri"
    year = 2001
    sites.append(site)
    locn = "Aberfeldy UK"
    df_leaf, g1 = get_leaf_data(df_wue, locn)
    df_leaf["site"] = site
    leaf_g1[site] = g1

    fname = os.path.join(fdir, "%s.%s.synth.hourly.allvars.csv" % (site, year))
    df_flux = pd.DataFrame()
    df, g1 = F.main(fname)
    df_flux = pd.concat([df_flux, df])
    flux_g1[site] = g1

    df_flux["gs_pred"] = FM.gs_model(df_flux["VPD_f"], df_flux["GPP_f"],
                                     df_flux["CO2"], 0.0, g1)
    df_flux["site"] = site

    df_leaf_all = pd.concat([df_leaf_all, df_leaf])
    df_flux_all = pd.concat([df_flux_all, df_flux])



    #############
    # Tumbarumba
    #############
    site = "AU-Tum"
    year = 2002
    sites.append(site)
    locn = "Tumbarumba flux tower Snowy Mts NSW"
    df_leaf, g1 = get_leaf_data(df_wue, locn)
    df_leaf["site"] = site
    leaf_g1[site] = g1

    fname = os.path.join(fdir, "%s.%s.synth.hourly.allvars.csv" % (site, year))
    df_flux = pd.DataFrame()
    df, g1 = F.main(fname)
    df_flux = pd.concat([df_flux, df])
    flux_g1[site] = g1

    df_flux["gs_pred"] = FM.gs_model(df_flux["VPD_f"], df_flux["GPP_f"],
                                     df_flux["CO2"], 0.0, g1)
    df_flux["site"] = site

    df_leaf_all = pd.concat([df_leaf_all, df_leaf])
    df_flux_all = pd.concat([df_flux_all, df_flux])



    #############
    # Puechabon
    #############
    site = "FR-Pue"
    year = 2006
    sites.append(site)
    locn = "Puechabon_France"
    df_leaf, g1 = get_leaf_data(df_wue, locn)
    df_leaf["site"] = site
    leaf_g1[site] = g1

    fname = os.path.join(fdir, "%s.%s.synth.hourly.allvars.csv" % (site, year))
    df_flux = pd.DataFrame()
    df, g1 = F.main(fname)
    df_flux = pd.concat([df_flux, df])
    flux_g1[site] = g1

    df_flux["gs_pred"] = FM.gs_model(df_flux["VPD_f"], df_flux["GPP_f"],
                                     df_flux["CO2"], 0.0, g1)
    df_flux["site"] = site

    df_leaf_all = pd.concat([df_leaf_all, df_leaf])
    df_flux_all = pd.concat([df_flux_all, df_flux])



    #############
    # Guyana
    #############
    site = "GF-Guy"
    year = 2006
    sites.append(site)
    locn = "Paracou_French Guiana"
    df_leaf, g1 = get_leaf_data(df_wue, locn)
    df_leaf["site"] = site
    leaf_g1[site] = g1

    fname = os.path.join(fdir, "%s.%s.synth.hourly.allvars.csv" % (site, year))
    df_flux = pd.DataFrame()
    df, g1 = F.main(fname)
    df_flux = pd.concat([df_flux, df])
    flux_g1[site] = g1

    df_flux["gs_pred"] = FM.gs_model(df_flux["VPD_f"], df_flux["GPP_f"],
                                     df_flux["CO2"], 0.0, g1)
    df_flux["site"] = site

    df_leaf_all = pd.concat([df_leaf_all, df_leaf])
    df_flux_all = pd.concat([df_flux_all, df_flux])


    #############
    # Harvard
    #############
    site = "US-Ha1"
    year = 1992
    sites.append(site)
    locn = "Prospect Hill Tract (Harvard Forest) USA"
    df_leaf, g1 = get_leaf_data(df_wue, locn)
    df_leaf["site"] = site
    leaf_g1[site] = g1

    fname = os.path.join(fdir, "%s.%s.synth.hourly.allvars.csv" % (site, year))
    df_flux = pd.DataFrame()
    df, g1 = F.main(fname)
    df_flux = pd.concat([df_flux, df])
    flux_g1[site] = g1

    df_flux["gs_pred"] = FM.gs_model(df_flux["VPD_f"], df_flux["GPP_f"],
                                     df_flux["CO2"], 0.0, g1)
    df_flux["site"] = site

    df_leaf_all = pd.concat([df_leaf_all, df_leaf])
    df_flux_all = pd.concat([df_flux_all, df_flux])




    #############
    # Soroe
    #############
    site = "DK-Sor"
    year = 1999
    sites.append(site)
    locn = "Soroe Denmark"
    df_leaf, g1 = get_leaf_data(df_wue, locn)
    df_leaf["site"] = site
    leaf_g1[site] = g1

    fname = os.path.join(fdir, "%s.%s.synth.hourly.allvars.csv" % (site, year))
    df_flux = pd.DataFrame()
    df, g1 = F.main(fname)
    df_flux = pd.concat([df_flux, df])
    flux_g1[site] = g1

    df_flux["gs_pred"] = FM.gs_model(df_flux["VPD_f"], df_flux["GPP_f"],
                                     df_flux["CO2"], 0.0, g1)
    df_flux["site"] = site

    df_leaf_all = pd.concat([df_leaf_all, df_leaf])
    df_flux_all = pd.concat([df_flux_all, df_flux])




    # Make plot
    colour_list = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


    fig = plt.figure(figsize=(12,9))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.1)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
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

    cnt = 0

    pfts = {'GF-Guy':'TRF', 'AU-Tum':'EBF', 'FR-Pue':'EBF', \
            'US-Ha1':'DBF', 'DK-Sor':'DBF', 'FI-Hyy':'ENF', \
            'FR-LBr':'ENF', 'UK-Gri':'ENF'}
    axes_index = [1,2,3,4,5,6,7,8]
    labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"]
    for site in sites:

        ax = fig.add_subplot(3,3,axes_index[cnt])

        leaf = df_leaf_all[df_leaf_all.site == site]
        flux = df_flux_all[df_flux_all.site == site]

        ax.plot(leaf.gs_index, leaf.gs, ls=" ", marker="o", markersize=2,
                 color=colour_list[0], markeredgecolor=colour_list[0],
                 alpha=0.6, zorder=1)
        p1, = ax.plot(leaf.gs_index, leaf.gs_pred, ls=" ", marker="o",
                 markersize=4, markeredgecolor="black", color=colour_list[0],
                 zorder=2, lw=4, label="Leaf gas exchange" )

        ax.plot(flux.gs_index, flux.gs, ls=" ",  marker="o", markersize=2,
                 color=colour_list[2], markeredgecolor=colour_list[1],
                 alpha=0.6, zorder=1)
        p2, = ax.plot(flux.gs_index, flux.gs_pred, ls=" ", marker="o",
                 markersize=4, markeredgecolor="black", color=colour_list[1],
                 zorder=2, lw=4, label="FLUXNET")



        ax.locator_params(nbins=4, axis="x")
        ax.locator_params(nbins=5, axis="y")
        ax.set_ylim(0, 2.5)


        plt.xticks([0.0, 0.05, 0.10, 0.15, 0.2], [0, 0.05, 0.1, 0.15, 0.2])
        ax.set_xlim(0, 0.2)

        if cnt != 0 and cnt != 3 and cnt != 6:
            plt.setp(ax.get_yticklabels(), visible=False)

        if cnt < 5:
            plt.setp(ax.get_xticklabels(), visible=False)


        if cnt == 7:
            ax.set_xlabel("Stomatal index (-)")
        if cnt == 3:
            ax.set_ylabel("$g_s$ (mol m$^{-2}$ s$^{-1}$)", position=(0.5,0.5))

        if cnt == 7:
            ax.legend([p1, p2], ["Leaf gas exchange", "FLUXNET"],
                      numpoints=1, loc="best", bbox_to_anchor=(1.4, 0.28))



        props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
        ax.text(0.03, 0.96, "%s %s (%s)" % (labels[cnt],site,pfts[site]),
                transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)

        ax.text(0.03, 0.85,
                "$g_{\mathrm{1,leaf}}$ = %.2f" % (leaf_g1[site]),
                transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)

        ax.text(0.03, 0.72,
                "$g_{\mathrm{1,flux}}$ = %.2f" % (flux_g1[site]),
                transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)


        cnt +=1



    odir = "/Users/%s/Dropbox/g1_leaf_ecosystem_paper/figures/figs/" % \
            (os.getlogin())
    #odir = "/Users/%s/Desktop/" % (os.getlogin())
    plt.savefig(os.path.join(odir, "g1_matching_leaf_and_canopy.pdf"),
                bbox_inches='tight', pad_inches=0.1)



class GetFluxData(FitFluxnetData): #subclass, inherits from FitFluxnetData

    def main(self, fname):
        df_site = pd.read_csv(self.site_fname)
        df_lai = pd.read_csv(self.lai_fname)

        d = self.get_site_info(df_site, df_lai, fname)
        df = pd.read_csv(fname, index_col='date',
                         parse_dates={'date': ["Year","DoY","Time"]},
                         date_parser=self.date_converter)

        # files contain a rouge date from the following year, fix it.
        df = self.fix_rogue_date(df, drop=True)
        (months, missing_gpp) = self.get_three_most_productive_months(df)

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

        # Filter extreme gs values which have come from extremely low VPD
        extreme = df['gs_est'].mean() + (3.0 * df['gs_est'].std())
        df = df[df['gs_est'] < extreme]

        F = FitMedlyn(fluxnet=True)
        params = F.setup_model_params()
        (result, success) = F.minimise_params(params, df, df["gs_est"])

        g1 = result.params['g1'].value
        g1_se = result.params['g1'].stderr
        g0 = 0.0
        df["gs_pred"] = F.gs_model(df["VPD_f"], df["GPP_f"], df["CO2"], g0, g1)

        df["gs_index"] = df['GPP_f'] / (df['CO2'] * np.sqrt(df['VPD_f']))

        # Drop the rest of the dataframe
        #df = df.rename(columns={'gs_est': 'gs'})
        df = df.rename(columns={'gs_est': 'gs'})

        return df, g1

def get_leaf_data(df_all, locn):

    df_out = pd.DataFrame()
    df_all = df_all.rename(columns={'Cond': 'gs'})

    L = FitMedlyn(fluxnet=False)

    by_PFT_list = ['Aobayama Sendai Japan','Paracou_French Guiana',\
                   'Tambopata_Peru']

    # Fit by: Location, then person, then species. Drop species if no. of
    # measurements < 6. unless location is one of above, then fit by
    # location

    df_locn = df_all[df_all.Location == locn]
    if df_locn.Location.iloc[0] in by_PFT_list:
        for pft in np.unique(df_locn['PFT']):
            df = df_locn[df_locn["PFT"] == pft]

            params = L.setup_model_params()
            (result, success) = L.minimise_params(params, df, df["gs"])
            g1 = result.params['g1'].value
            g1_se = result.params['g1'].stderr
            g0 = 0.0

            df_out = pd.concat([df_out,df])
            df_out["gs_pred"] = L.gs_model(df_out["VPD"], df_out["Photo"],
                                           df_out["CO2S"], g0, g1)
    else:
        g1_sum = 0.0
        n = 0
        for ppl in np.unique(df_locn['Datacontrib']):
            df_ppl = df_locn[df_locn.Datacontrib == ppl]
            for spp in np.unique(df_ppl['Species']):
                df = df_ppl[df_ppl.Species == spp]
                if (len(df["gs"]) > 5):
                    params = L.setup_model_params()
                    (result, success) = L.minimise_params(params, df, df["gs"])
                    g1 = result.params['g1'].value
                    g1_se = result.params['g1'].stderr
                    g0 = 0.0
                    g1_sum += g1
                    n += 1
                    df_out = pd.concat([df_out,df])
        g1 = g1_sum / float(n)
        df_out["gs_pred"] = L.gs_model(df_out["VPD"], df_out["Photo"],
                                       df_out["CO2S"], g0, g1)


    df_out["gs_index"] = df_out['Photo'] / (df_out['CO2S'] * \
                            np.sqrt(df_out['VPD']))

    return df_out, g1

def add_pft_info(df):
    """
    Add missing PFT categories to dataframe
    """
    df.loc[:,'PFT'] = pd.Series(["missing"]*len(df), index=df.index)
    df.loc[:,'Scale'] = pd.Series(["Leaf gas exchange"]*len(df),
                                  index=df.index)
    for index, row in df.iterrows():
        if (row.Type == 'gymnosperm' and
            row.Plantform == 'tree' and
            row.Leafspan == 'evergreen'):
            df.loc[index,'PFT'] = "ENF"
        elif (row.Type == 'angiosperm' and
              row.Plantform == 'tree' and
              row.Leafspan == 'deciduous'):
              df.loc[index,'PFT'] = "DBF"
        elif (row.Type == 'angiosperm' and
              row.Leafspan == 'evergreen' and
              row.Plantform == 'tree' and
              row.Tregion == 'temperate'):
              df.loc[index,'PFT'] = "EBF"
        elif (row.Type == 'angiosperm' and
              row.Leafspan == 'evergreen' and
              row.Plantform == 'tree' and
              row.Tregion == 'tropical'):
              df.loc[index,'PFT'] = "TRF"
        elif (row.Type == 'angiosperm' and
              row.Plantform == 'shrub'):
              df.loc[index,'PFT'] = "SHB"
        elif (row.Pathway == 'C3' and
              row.Plantform == 'crop'):
              df.loc[index,'PFT'] = "CRO"
        elif (row.Pathway == 'C3' and
              row.Plantform == 'grass'):
              df.loc[index,'PFT'] = "C3G"
        elif (row.Pathway == 'C4' and
              row.Plantform == 'grass'):
              df.loc[index,'PFT'] = "C4G"
        elif (row.Plantform == 'savanna' and
              row.Leafspan == 'evergreen'):
              df.loc[index,'PFT'] = "SAV"
        elif (row.Plantform == 'savanna' and
              row.Leafspan == 'deciduous'):
              df.loc[index,'PFT'] = "SAV"

    return (df)

if __name__ == "__main__":

    main()
