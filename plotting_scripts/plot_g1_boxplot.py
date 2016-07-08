#!/usr/bin/env python

"""
Plot Fluxnet fits by PFT

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (05.11.2015)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import brewer2mpl

def main():

    fdir = "data/processed"
    df_flux = pd.read_csv(os.path.join(fdir, "g1_fluxnet_screened.csv"))
    df_leaf = pd.read_csv(os.path.join(fdir, "g1_leaf_gas_exchange.csv"))
    df_isotope = pd.read_csv(os.path.join(fdir, "g1_isotope_screened.csv"))
    df_sig = pd.read_csv(os.path.join(fdir, "g1_siglet_withinPFT_coupled.csv"))

    df = pd.concat([df_flux, df_leaf, df_isotope])

    # Estimate data loss due to cut off...
    """
    aa = len(df_flux)
    df_flux = df_flux[df_flux.g1 <= 14]
    bb = len(df_flux)
    print (aa - bb) / float(aa) * 100.

    aa = len(df_leaf)
    df_leaf = df_leaf[df_leaf.g1 <= 14]
    bb = len(df_leaf)
    print (aa - bb) / float(aa) * 100.

    aa = len(df_isotope)
    df_isotope = df_isotope[df_isotope.g1 <= 14]
    bb = len(df_isotope)
    print (aa - bb) / float(aa) * 100.
    """
    # Cut off ylimit at 14 for visual purposes
    aa = len(df)
    df = df[df.g1 <= 14]
    bb = len(df)
    #print (aa - bb) / float(aa) * 100.

    pfts = ['ENF','EBF','DBF','TRF','SAV','SHB','C3G','C4G','C3C','C4C']

    sns.set_style("ticks")
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


    fig = plt.figure(figsize=(12,6))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
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

    ax = fig.add_subplot(111)

    # define outlier properties
    flierprops = dict(marker='o', markersize=3, markerfacecolor="black")
    ax = sns.boxplot(x="PFT", y="g1", hue="Scale", data=df, palette="Set2",
                     order=pfts, flierprops=flierprops,
                     hue_order=['Leaf gas exchange','Leaf isotope','FLUXNET'])

    ax.set_ylabel("$g_1$ (kPa$^{0.5}$)")
    ax.set_xlabel("Plant functional type")

    ax.set_ylim(-1.5, 14)
    ax.set_yticks(np.arange(0,16,2))


    ax.set_xticklabels(pfts)

    for i,p in enumerate(pfts[:-3]):

        if i == 0:
            offset = 0

        sig = df_sig[(df_sig.PFT  == p) &
                     (df_sig.Method == "leaf")].siglet.values[0]

        if p == "EBF":
            plt.text(-0.31+offset-0.04, -1,"%s\n" % (sig),
                     horizontalalignment='left', size='x-small', color="k")
        else:
            plt.text(-0.31+offset, -1,"%s\n" % (sig),
                     horizontalalignment='left', size='x-small', color="k")
        sig = df_sig[(df_sig.PFT  == p) &
                     (df_sig.Method == "isotope")].siglet.values[0]
        plt.text(-0.04+offset, -1,"%s\n" % (sig), horizontalalignment='left',
                 size='x-small', color="k")
        sig = df_sig[(df_sig.PFT == p) &
                     (df_sig.Method == "fluxnet")].siglet.values[0]
        plt.text(0.23+offset, -1,"%s\n" % (sig), horizontalalignment='left',
                 size='x-small', color="k")

        offset += 1.0

    for i,p in enumerate(pfts):

        if i == 0:
            offset = 0

        if p!= "C4C":
            plt.text(-0.13+offset, 15+4,
                     "n=%d\n" % (len(df_leaf[df_leaf.PFT == p].g1)),
                     horizontalalignment='left', size='small', color=set2[0],
                     weight="bold")

        if p != "C4G" and p != "C3C" and p!= "C4C":
            plt.text(-0.13+offset, 14.5+4,
                     "n=%d\n" % (len(df_isotope[df_isotope.PFT == p].g1)),
                     horizontalalignment='left', size='small', color=set2[1],
                     weight="bold")

        if p != "C4G":
            plt.text(-0.13+offset, 14+4,
                     "n=%d\n" % (len(df_flux[df_flux.PFT  == p].g1)),
                     horizontalalignment='left', size='small', color=set2[2],
                     weight="bold")

        offset += 1.0

    plt.legend(loc="upper right")

    odir = "/Users/%s/Dropbox/g1_leaf_ecosystem_paper/figures/figs/" % \
            (os.getlogin())
    plt.savefig(os.path.join(odir,
                "g1_gas_exchange_isotope_fluxnet_boxplot.pdf"),
                bbox_inches='tight', pad_inches=0.1)

    #plt.show()

if __name__ == "__main__":

    main()
