#!/usr/bin/env python

"""
Plot g1 vs. C4 frac

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (18.04.2016)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import brewer2mpl

def main():

    fdir = "data/processed"
    df_flux = pd.read_csv(os.path.join(fdir, "g1_fluxnet_screened.csv"))

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

    colour_list = brewer2mpl.get_map('Accent', 'qualitative', 8).mpl_colors
    # CB palette  with grey:
    # from http://jfly.iam.u-tokyo.ac.jp/color/image/pallete.jpg
    #colour_list = ["#56B4E9", "#009E73", "#0072B2", "#F0E442",\
    #               "#E69F00", "#D55E00", "#CC79A7", "#999999"]


    ax1 = fig.add_subplot(111)

    flux = df_flux[(df_flux.PFT == "C3G")]
    flux = flux[flux.g1 <= 18] # cut off to match fig. 1

    ax1.plot(flux.c4_frac, flux.g1,  ls=" ", marker="o", color="black",
             markeredgecolor="lightgrey", alpha=0.8)

    ax1.set_xlim(-0.05, 0.45)

    ax1.locator_params(nbins=7, axis="x")

    ax1.set_xlabel("C4 fraction (-)")
    ax1.set_ylabel("Estimated $g_1$ (kPa$^{0.5}$)")


    odir = "/Users/%s/Dropbox/g1_leaf_ecosystem_paper/figures/figs/" % \
            (os.getlogin())
    plt.savefig(os.path.join(odir, "g1_vs_c4frac.pdf"),
                bbox_inches='tight', pad_inches=0.1)

    #plt.show()

if __name__ == "__main__":

    main()
