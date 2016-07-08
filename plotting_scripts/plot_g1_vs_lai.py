#!/usr/bin/env python

"""
Plot g1 fluxnet vs. LAI for forest and non-forest vegetation

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (17.11.2015)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import brewer2mpl
import seaborn as sns

def main():

    fdir = "data/processed"
    df_flux = pd.read_csv(os.path.join(fdir, "g1_fluxnet_screened.csv"))
    df = df_flux[(df_flux.lai_max > 0.0)]
    df = df[(df.lai_max < 100.0)]

    sns.set_style("ticks")
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    golden_mean = 0.6180339887498949
    width = 9
    height = width * golden_mean
    fig = plt.figure(figsize=(width, height))
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

    ax1 = fig.add_subplot(111)
    forest_PFTs = ['ENF', 'EBF', 'DBF', 'TRF']
    nonforest_PFTs = ['SAV', 'SHB', 'C3G']

    forest = df[df.PFT.isin(forest_PFTs)]
    nonforest = df[df.PFT.isin(nonforest_PFTs)]


    ax1.plot(forest.lai_max,forest.g1, ls=" ", marker="o", color="darkgreen",
             label='Forest', alpha=0.9)

    ax1.plot(nonforest.lai_max,nonforest.g1, ls=" ", marker="o",
             color="lightgreen", label='Non-Forest', alpha=0.9)

    ax1.legend(numpoints=1, ncol=1, loc="best", frameon=False)


    ax1.locator_params(nbins=7, axis="x")

    ax1.set_xlabel("LAI (m$^{2}$ m$^{-2}$)")
    ax1.set_ylabel("Estimated $g_1$ (kPa$^{0.5}$)")


    odir = "/Users/%s/Dropbox/g1_leaf_ecosystem_paper/figures/figs/" % \
            (os.getlogin())
    plt.savefig(os.path.join(odir, "g1_vs_lai.pdf"),
                    bbox_inches='tight', pad_inches=0.1)

    #plt.show()

if __name__ == "__main__":

    main()
