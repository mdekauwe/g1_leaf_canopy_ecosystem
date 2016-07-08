#!/usr/bin/env python

"""
Put FLUXNET, gas exchange and isotope site locations onto a map

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (29.04.2016)"
__email__ = "mdekauwe@gmail.com"


from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import glob

import brewer2mpl
import numpy as np
def main():

    fdir = "data/processed"
    df_flux = pd.read_csv(os.path.join(fdir, "g1_fluxnet_screened.csv"))
    df_leaf = pd.read_csv(os.path.join(fdir, "g1_leaf_gas_exchange.csv"))
    df_isotope = pd.read_csv(os.path.join(fdir, "g1_isotope_screened.csv"))

    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
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


    set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors

    golden_mean = 0.6180339887498949
    width = 7
    height = width * golden_mean
    fig = plt.figure(figsize=(width, height))
    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.00)
    #m = Basemap(resolution='c', projection='robin', lon_0=0)
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
    m.fillcontinents(color='lightgray',lake_color='white',zorder=1)
    m.drawcoastlines(linewidth=0.05, color='k')
    m.drawcountries(linewidth=0.1, color='k')
    m.drawparallels(np.array([-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0]),
                    labels=[1,0,0,0], fontsize=10, linewidth=0.0)
    m.drawmeridians(np.array([-180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0]),
                    labels=[0,0,0,1], fontsize=10, linewidth=0.0)

    # compute the native map projection coordinates for locations.
    x1,y1 = m(df_flux.longitude.values, df_flux.latitude.values)
    x2,y2 = m(df_isotope.longitude.values, df_isotope.latitude.values)
    x3,y3 = m(df_leaf.longitude.values, df_leaf.latitude.values)


    m.scatter(x1, y1, color=set2[2], s=7.5, label="FLUXNET", zorder=2,
              alpha=0.7)
    m.scatter(x2, y2, color=set2[1], s=7.5, label="Leaf isotope", zorder=2,
              alpha=0.7)
    m.scatter(x3, y3, color=set2[0], s=7.5, label="Leaf gas exchange", zorder=2,
              alpha=0.7)

    lgnd = plt.legend(scatterpoints=1, loc=(0.0, 0.15), frameon=False)

    # Increase the size of legend markers
    lgnd.legendHandles[0]._sizes = [40]
    lgnd.legendHandles[1]._sizes = [40]
    lgnd.legendHandles[2]._sizes = [40]

    #m.bluemarble()
    #m.shadedrelief()
    #plt.show()
    odir = "/Users/%s/Dropbox/g1_leaf_ecosystem_paper/figures/figs/" % \
            (os.getlogin())
    #odir = "/Users/%s/Desktop/" % (os.getlogin())
    plt.savefig(os.path.join(odir,
                "location_map.png"),
                bbox_inches='tight', pad_inches=0.1, dpi=600)


if __name__ == "__main__":

    main()
