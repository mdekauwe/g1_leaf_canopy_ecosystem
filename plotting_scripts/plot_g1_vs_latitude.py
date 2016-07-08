#!/usr/bin/env python

"""
Plot g1 leaf gas exchange, fluxnet, isotope vs. latitude.

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
    df_leaf = pd.read_csv(os.path.join(fdir, "g1_leaf_gas_exchange.csv"))
    df_isotope = pd.read_csv(os.path.join(fdir, "g1_isotope_screened.csv"))

    sns.set_style("ticks")
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    fig = plt.figure(figsize=(12,12))
    fig.subplots_adjust(hspace=0.05)
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

    colour_list = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors

    #colour_list = sns.palplot(sns.color_palette("colorblind", 10))
    #colour_list = sns.color_palette("Set2", 10)
    # CB palette  with grey:
    # from http://jfly.iam.u-tokyo.ac.jp/color/image/pallete.jpg
    #colour_list = ["#56B4E9", "#009E73", "#0072B2", "#F0E442",\
    #               "#E69F00", "#D55E00", "#CC79A7", "#999999"]
    colour_list = sns.color_palette("Accent", 10)

    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)

    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    pft_order = ['ENF', 'EBF', 'DBF', 'TRF']
    for i, pft in enumerate(pft_order):
        leaf = df_leaf[df_leaf.PFT == pft]
        isotope = df_isotope[df_isotope.PFT == pft]
        flux = df_flux[df_flux.PFT == pft]


        if i >= 3:
            cidx = i+1
        else:
            cidx = i

        for lat in np.unique(leaf.latitude):
            data_lat = leaf[leaf.latitude == lat]

            ax1.errorbar(np.mean(data_lat.g1), np.mean(data_lat.latitude),
                         xerr=stats.sem(data_lat.g1), ls=" ", marker="o",
                         color=colour_list[cidx], markeredgecolor="lightgrey",
                         alpha=0.8, capsize=False)

        for lat in np.unique(isotope.latitude):
            data_lat = isotope[isotope.latitude == lat]

            ax2.errorbar(np.mean(data_lat.g1), np.mean(data_lat.latitude),
                         xerr=stats.sem(data_lat.g1), ls=" ", marker="o",
                         color=colour_list[cidx], markeredgecolor="lightgrey",
                         alpha=0.8, capsize=False)

        for lat in np.unique(flux.latitude):
            data_lat = flux[flux.latitude == lat]

            ax3.errorbar(np.mean(data_lat.g1), np.mean(data_lat.latitude),
                         xerr=stats.sem(data_lat.g1), ls=" ", marker="o",
                         color=colour_list[cidx], markeredgecolor="lightgrey",
                         alpha=0.8, capsize=False)

    for i, pft in enumerate(pft_order):
        if i >= 3:
            cidx = i+1
        else:
            cidx = i

        ax1.plot(np.nan, np.nan, ls=" ", marker="o", color=colour_list[cidx],
                 markeredgecolor="lightgrey", label=pft, alpha=0.9)
    ax1.legend(numpoints=1, ncol=1, loc="best", frameon=False)

    pft_order = ['SAV', 'SHB', 'C3G', 'C4G', 'C3C', 'C4C']
    for i, pft in enumerate(pft_order):
        leaf = df_leaf[df_leaf.PFT == pft]
        isotope = df_isotope[df_isotope.PFT == pft]
        flux = df_flux[df_flux.PFT == pft]

        if i >= 3:
            cidx = i+1
        else:
            cidx = i
        for lat in np.unique(leaf.latitude):
            data_lat = leaf[leaf.latitude == lat]

            ax4.errorbar(np.mean(data_lat.g1), np.mean(data_lat.latitude),
                         xerr=stats.sem(data_lat.g1), ls=" ", marker="D",
                         color=colour_list[cidx], markeredgecolor="lightgrey",
                         alpha=0.8, capsize=False)


        for lat in np.unique(isotope.latitude):
            data_lat = isotope[isotope.latitude == lat]

            ax5.errorbar(np.mean(data_lat.g1), np.mean(data_lat.latitude),
                         xerr=stats.sem(data_lat.g1), ls=" ", marker="D",
                         color=colour_list[cidx], markeredgecolor="lightgrey",
                         label=pft, alpha=0.8, capsize=False)
        for lat in np.unique(flux.latitude):
            data_lat = flux[flux.latitude == lat]

            ax6.errorbar(np.mean(data_lat.g1), np.mean(data_lat.latitude),
                         xerr=stats.sem(data_lat.g1), ls=" ", marker="D",
                         color=colour_list[cidx], markeredgecolor="lightgrey",
                         alpha=0.8, capsize=False)

    for i, pft in enumerate(pft_order):
        if i >= 3:
            cidx = i+1
        else:
            cidx = i

        ax4.plot(np.nan, np.nan, ls=" ", marker="D", color=colour_list[cidx],
                 markeredgecolor="lightgrey", label=pft, alpha=0.9)


    ax4.legend(numpoints=1, ncol=1, loc="best", frameon=False)

    labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")


    for i, ax in enumerate([ax1, ax2, ax3, ax4, ax5, ax6]):
        ax.set_xlim(0, 14)
        ax.set_ylim(-60, 90)
        ax.locator_params(nbins=6, axis="x")
        ax.locator_params(nbins=6, axis="y")

        ax.axhline(y=0.0, c='grey', lw=1.0, ls='--')
        ax.axhline(y=-23.43723, c='grey', lw=1.0, ls='-.')
        ax.axhline(y=23.43723, c='grey', lw=1.0, ls='-.')

        ax.text(0.03, 0.98, labels[i], transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)

    for ax in [ax1, ax2, ax3]:
        plt.setp(ax.get_xticklabels(), visible=False)

    for ax in [ax2, ax3, ax5, ax6]:
        plt.setp(ax.get_yticklabels(), visible=False)

    ax1.set_title("Leaf gas exchange")
    ax2.set_title("Leaf isotope")
    ax3.set_title("FLUXNET")


    ax1.set_ylabel("Latitude (degrees)", position=(0.5, 0.0))
    ax5.set_xlabel("Estimated $g_1$ (kPa$^{0.5}$)")


    odir = "/Users/%s/Dropbox/g1_leaf_ecosystem_paper/figures/figs/" % \
            (os.getlogin())
    plt.savefig(os.path.join(odir, "g1_vs_latitude.pdf"),
                bbox_inches='tight', pad_inches=0.1)

    #plt.show()

if __name__ == "__main__":

    main()
