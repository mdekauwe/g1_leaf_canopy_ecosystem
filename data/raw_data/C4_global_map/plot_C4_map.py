#!/usr/bin/env python

"""
Plot NACP Synergetic land cover product (SYNMAP) Global C3 and C4 maps
0.5-degree: present relative fraction of C4 grasses v2.0

Reference:
----------
Jung, M., Henkel, K., Herold, M., and Churkina, G.: Exploiting synergies of
global land cover products for carbon cycle modeling, Remote Sens. Environ.,
101, 534–553, doi:10.1016/j.rse.2006.01.020, 2006.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (15.04.2016)"
__email__ = "mdekauwe@gmail.com"

import numpy as np
import matplotlib.pyplot as plt
import gdal
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import AxesGrid
import glob
import os
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sns
import pandas as pd

def colorbar_index(cax=None, ncolours=None, cmap=None, orientation=None):
    cmap = cmap_discretize(cmap, ncolours)
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolours+0.5)
    colorbar = plt.colorbar(mappable, cax=cax, orientation=orientation)
    colorbar.set_ticks(np.linspace(0, ncolours, ncolours))

    return colorbar

def cmap_discretize(cmap, N):
    """
    Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)

    """
    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
                      for i in xrange(N + 1)]
    return mcolors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)

# load fluxnet sites
fname = '../../processed/g1_fluxnet_screened.csv'
df = pd.read_csv(fname, sep=",")

x = []
y = []
for i, row in df.iterrows():
    y.append(row.latitude)
    x.append(row.longitude)

# open dataset
ds = gdal.Open('wcs.tiff')
data = np.array(ds.GetRasterBand(1).ReadAsArray())
ds = None

fig = plt.figure(figsize=(14, 10))
grid = AxesGrid(fig, [0.05,0.05,0.9,0.9], nrows_ncols=(1,1), axes_pad=0.1,
                cbar_mode='single', cbar_pad=0.2, cbar_size="3%",
                cbar_location='bottom', share_all=True)

m = Basemap(projection='cyl', llcrnrlon=-180.0, llcrnrlat=-90.0, \
            urcrnrlon=180, urcrnrlat=90.0, resolution='c')



# Range on colourbar
ncolours = 11
vmin = 0.0
vmax = 1.0

bmap = sns.blend_palette(["white", "darkgreen"], ncolours, as_cmap=True)
ax = grid[0]
m.ax = ax
m.drawcoastlines(linewidth=0.1, color='k')
m.drawcountries(linewidth=0.1, color='k')
image = m.imshow(np.flipud(data), bmap,
                 colors.Normalize(vmin=vmin, vmax=vmax, clip=True),
                 interpolation='nearest')
cbar = colorbar_index(cax=grid.cbar_axes[0], ncolours=ncolours, cmap=bmap,
                      orientation='horizontal')
cbar.set_ticklabels(np.linspace(vmin, vmax, ncolours))
cbar.set_label("C4 Fraction (-)", fontsize=16)

# fluxnet sites
x1, y1 = m(x, y)
m.scatter(x1, y1, marker="o", color="black", alpha=0.7)

fig.savefig("C4_frac.png", bbox_inches='tight', pad_inches=0.1, dpi=300)
plt.show()

# close dataset
