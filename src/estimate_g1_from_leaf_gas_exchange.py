#!/usr/bin/env python

"""
Fit Belinda's stomatal conductance model to gas exchange data to obtain
estimates of g1

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (16.04.2016)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from rmse import rmse

from fit_medlyn_gs_model import FitMedlyn

class FitGasExchange(object):
    """
    Fit Belinda's stomatal conductance model using the lmfit package to leaf
    gas exchange data
    """

    def __init__(self, fdir, ofdir, fname, ofname):

        self.fname = os.path.join(fdir, fname)
        self.ofname = os.path.join(ofdir, ofname)
        self.out_cols = ["Scale","g1","g1_se","r2","RMSE","n","Species",\
                         "Datacontrib","Location","latitude","longitude","PFT",\
                         "Pathway","Type","Plantform", "Leafspan","Tregion"]

    def main(self):
        df_all = pd.read_csv(self.fname)
        df_all = self.add_pft_info(df_all)

        # set up output file
        df_out = pd.DataFrame(columns=self.out_cols)

        # Omit Phragmites - it is a wetland plant
        df_all = df_all[df_all.Species != "Phragmites communis"]

        # List of locations where there is insufficient data by species so fit
        # by PFT
        by_PFT_list = ['Aobayama Sendai Japan','Paracou_French Guiana',\
                       'Tambopata_Peru']

        # Fit by: Location, then person, then species. Drop species if no. of
        # measurements < 6. unless location is one of above, then fit by
        # location
        for locn in np.unique(df_all['Location']):
            df_locn = df_all[df_all.Location == locn]
            if df_locn.Location.iloc[0] in by_PFT_list:
                for pft in np.unique(df_locn['PFT']):
                    df = df_locn[df_locn["PFT"] == pft]
                    df_out = self.fitting_wrapper(df, df_out)
            else:
                for ppl in np.unique(df_locn['Datacontrib']):
                    df_ppl = df_locn[df_locn.Datacontrib == ppl]
                    for spp in np.unique(df_ppl['Species']):
                        df = df_ppl[df_ppl.Species == spp]
                        if (len(df["Cond"]) > 5):
                            df_out = self.fitting_wrapper(df, df_out)

        # drop any nan data
        df_out = df_out.drop(df_out[np.isinf(df_out.g1_se)].index)

        # screen bad fits in a consistent way with the fluxnet data
        df_out = df_out.drop(df_out[df_out.r2 < 0.2].index)

        df_out.to_csv(self.ofname, index=False)

    def fitting_wrapper(self, df, df_out):
        F = FitMedlyn(fluxnet=False)
        params = F.setup_model_params()
        (result, success) = F.minimise_params(params, df, df["Cond"])

        if success:
            g1 = result.params['g1'].value
            g1_se = result.params['g1'].stderr
            g0 = 0.0
            model = F.gs_model(df["VPD"], df["Photo"], df["CO2S"], g0, g1)
            rsq = (pearsonr(df["Cond"], model)[0])**2
            num_pts = len(df["Cond"])
            rmse_val = rmse(df["Cond"], model)

            row = pd.Series([df.Scale.iloc[0], g1, g1_se, rsq, rmse_val,
                             num_pts, df.Species.iloc[0],
                             df.Datacontrib.iloc[0], df.Location.iloc[0],
                             df.latitude.iloc[0], df.longitude.iloc[0],
                             df.PFT.iloc[0], df.Pathway.iloc[0],
                             df.Type.iloc[0], df.Plantform.iloc[0],
                             df.Leafspan.iloc[0], df.Tregion.iloc[0]],
                             index=self.out_cols)

            df_out = df_out.append(row, ignore_index=True)

        return (df_out)

    def add_pft_info(self, df):
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
                  df.loc[index,'PFT'] = "C3C"
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

    G = FitGasExchange(fdir="data/raw_data/leafgasexchange/",
                       ofdir="data/processed/",
                       fname="WUEdatabase_13_04_2016_fixed.csv",
                       ofname="g1_leaf_gas_exchange.csv")
    G.main()
