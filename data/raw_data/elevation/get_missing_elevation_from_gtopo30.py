# -*- coding: UTF-8 -*-
"""
Grabbing the missing elevation data from the GTOPO30 is a global digital
elevation model (DEM) with a horizontal grid spacing of 30 arc seconds
(approximately 1 kilometer).

http://www.geonames.org/export/ws-overview.html

Using a bit of JSON logic to query the API of the site.

"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (22.02.2016)"
__email__   = "mdekauwe@gmail.com"

import pandas as pd
import numpy as np
import sys
import requests

df_gap = pd.read_csv("gappy_elevation.csv")


new_elev = []
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

    if lat < -500.0 or lon < -500.0:
        new_elev.append(-999.9)
    elif df_gap["Elevation"][i] < -500.0:
        r = requests.get('http://api.geonames.org/gtopo30JSON?lat=%f&lng=%f&username=mdekauwe' % (lat, lon))
        #print r.json()
        elev = r.json()['gtopo30']
        #print lat, lon, elev
        new_elev.append(elev)
    else:
        new_elev.append(df_gap["Elevation"][i])


print "Missing elevation data gap-filled from gtopo30 1km data"
print "using get_missing_elevation_from_gtopo30.py"
print "Site_ID,Latitude,Longitude,Elevation"
for i in xrange(len(df_gap)):
    print "%s,%s,%s,%f" % (df_gap["Site_ID"][i], df_gap["Latitude"][i],\
                           df_gap["Longitude"][i],new_elev[i])
