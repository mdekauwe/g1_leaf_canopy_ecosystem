# -*- coding: UTF-8 -*-
"""
Grabbing the missing elevation data from the GTOPO30 is a global digital
elevation model (DEM) with a horizontal grid spacing of 30 arc seconds
(approximately 1 kilometer).

http://www.geonames.org/export/ws-overview.html

Using a bit of JSON logic to query the API of the site.

- Check our grabbed values match the ones we already have!

"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (22.02.2016)"
__email__   = "mdekauwe@gmail.com"

import matplotlib
matplotlib.use('agg') # stop windows popping up

import pandas as pd
import numpy as np
import sys
import requests
import matplotlib.pyplot as plt

df_gap = pd.read_csv("gappy_elevation.csv")

new_elev = []
orig_elev = []
for i in xrange(len(df_gap)):

    lat = df_gap["Latitude"][i].replace(" ", "")
    lon = df_gap["Longitude"][i].replace(" ", "")

    if lat == "TBD":
        lat = -999.9
    else:
        lat = float(lat)

    if lon == "TBD":
        lon = -999.9
    else:
        lon = float(lon)

    if lat > -500.0 and lon > -500 and df_gap["Elevation"][i] > -500:
        r = requests.get('http://api.geonames.org/gtopo30JSON?lat=%f&lng=%f&username=mdekauwe' % (lat, lon))
        #print r.json()
        elev = r.json()['gtopo30']
        #print lat, lon, df_gap["Elevation"][i], elev
        new_elev.append(elev)
        orig_elev.append(df_gap["Elevation"][i])
print "done"

fig = plt.figure(figsize=(6,6))
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

ax1 = fig.add_subplot(111, aspect='equal')
ax1.plot(orig_elev, new_elev, "ko", alpha=0.3)

ax1.set_ylabel("New (m)")
ax1.set_xlabel("Old (m)")
plt.savefig("elevation_comparison.pdf", bbox_inches='tight', pad_inches=0.1)
#plt.show()
