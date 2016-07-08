import pandas as pd
import sys
import numpy as np

years = ["1991", "1992", "1993", "1994", "1995", "1996", "1997", "1998", \
         "1999", "2000", "2001", "2002", "2003", "2004", "2005", \
         "2006", "2007"]
df = pd.read_csv("sites_use.csv")

n = 0
nn = 0

f = open("free_sites.csv", "w")
print >> f, "site,yr"
for i in xrange(len(df)):
    site = df["site"][i]
    for j, yr in enumerate(years):

        #if df[yr][i] == "Fair_Use" or df[yr][i] == "Open":
        if df[yr][i] == "Fair_Use":
            print >> f, "%s,%d" % (site, int(yr))

            n+=1
        elif df[yr][i] == "LaThuile_2007":
            nn+=1

f.close()

print n, nn
