#!/usr/bin/env python

"""
Estimate WUE from leaf C-isotope data. This approach is based on the the fact
that during photosynthesis, photosynthetic enzymes discriminate against the
heavier stable isotope 13C, relative to 12C. Consequently, C in leaves is always
depleted in 13C when compared to the atmosphere.

The extent of the discrimination depends on Ci. Where Ci is low relative to Ca,
the air inside the leaf is enriched with 13C and the ability of the enzyme to
discriminate declines. In this scenario plants fix a greater proportion of 13C
than plants performing photosynthesis at a higher Ci.


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
from lmfit import minimize, Parameters
import matplotlib.pyplot as plt
import calendar
import datetime as dt
from scipy.stats import pearsonr
from rmse import rmse


class FitLeafIsotope(object):

    def __init__(self, fdir, ofdir, fname, pft_fname, vpd_fname, ofname):

        self.fname = os.path.join(fdir, fname)
        self.pft_fname = os.path.join(fdir, pft_fname)
        self.vpd_fname = os.path.join(fdir, vpd_fname)
        self.ofname = os.path.join(ofdir, ofname)

    def main(self):

        df_isotope = pd.read_csv(self.fname)

        # PFT info comes from YS and Wang Han, they made some mistakes which
        # I have fixed where possible
        df_pft = pd.read_csv(self.pft_fname)

        df_vpd = pd.read_csv(self.vpd_fname)

        df_isotope = self.add_missing_PFT_data(df_isotope, df_pft)
        df_isotope = self.add_missing_VPD_mGDD0_data(df_isotope, df_vpd)
        df_isotope = self.estimate_g1(df_isotope)

        df_isotope = df_isotope.drop(df_isotope[df_isotope['g1'] < 0.0].index)

        # Add tropical rain forest identifier
        idx = ( ((df_isotope.latitude >= -23.43723) &
                 (df_isotope.latitude <= 23.43723) &
                 (df_isotope.PFT == "EBF")) |
                ((df_isotope.latitude >= -23.43723) &
                 (df_isotope.latitude <= 23.43723) &
                 (df_isotope.PFT == "DBF")) |
                ((df_isotope.latitude >= -23.43723) &
                 (df_isotope.latitude <= 23.43723) &
                 (df_isotope.PFT == "ENF")) )
        df_isotope.loc[idx, "PFT"] = "TRF"

        # some reclassificaiton of labels
        df_isotope.loc[df_isotope.PFT == "CRO", "PFT"] = "C3C"

        df_isotope = df_isotope[df_isotope.PFT != "-9999"]

        # drop C4 stuff, or anything claiming to be c4
        df_isotope = df_isotope[df_isotope.PFT != "C4G"]

        df_isotope['Scale'] = 'Leaf isotope'

        df_isotope.to_csv(self.ofname, index=False)

    def add_missing_VPD_mGDD0_data(self, df_isotope, df_vpd):

        matching_vpd = []
        matching_mGDD0 = []
        for i in xrange(len(df_isotope)):

            # find closest match
            distance = (np.abs(df_vpd["lat"] - df_isotope["latitude"][i])) + \
                       (np.abs(df_vpd["lon"] - df_isotope["longitude"][i]))
            index = np.argmin(distance)

            matching_vpd.append(df_vpd["mVPD0_daytime"][index])
            matching_mGDD0.append(df_vpd["mGDD0"][index])

        df_isotope["VPD"] = matching_vpd
        df_isotope["mGDD0"] = matching_mGDD0

        return df_isotope

    def add_missing_PFT_data(self, df_isotope, df_pft):

        matching_pft = []
        for i in xrange(len(df_isotope)):

            study_ID = df_isotope["study_ID"][i]
            site_ID = df_isotope["site_ID"][i]
            orig_spp = df_isotope["orig_spp"][i]
            family = df_isotope["family"][i]
            lat = df_isotope["latitude"][i]
            lon = df_isotope["longitude"][i]

            found_pft = df_pft[(df_pft["study_ID"] == study_ID) &
                               (df_pft["site_ID"] == site_ID) &
                               (df_pft["orig_spp"] == orig_spp) &
                               (df_pft["family"] == family) &
                               (df_pft["latitude"] == lat) &
                               (df_pft["longitude"] == lon)].PFT.values[0]
            matching_pft.append(found_pft)
        df_isotope["PFT"] = matching_pft

        return df_isotope

    def estimate_g1(self, df):

        df["ci_ca"] = np.ones(len(df)) * -9999.9
        df["g1"] = np.ones(len(df)) * -9999.9

        for i in xrange(len(df)):

            if df["ps_type"][i] == "C3":
                df.loc[i,"ci_ca"] = self.delta_c3(df["big_delta_merged"][i])
            #elif df["ps_type"][i] == "C4":
            #    df.loc[i, "ci_ca"] = self.delta_c4(df["big_delta_merged"][i],
            #                                       df["mGDD0"][i])

            # missing is -9999
            if df["VPD"][i] >= 0.0:
                df.loc[i,"g1"] = self.calculate_g1(df["ci_ca"][i], df["VPD"][i])

        return df

    def delta_c3(self, delta):
        """
        delta (photosyntheric 13C discrimination) = a + (b - a) * Ci/Ca
        delta allows a time-integrated estimate of Ci:Ca and intrinsic WUE.

        References:
        -----------
        * Farquhar GD, O'Leary MH, Berry JA (1982) On the relationship
          between carbon isotope discrimination and the intercellular
          carbon dioxide concentration in leaves. Aust J Plant Physiol
          9:121-137
        """
        a = 4.4   # fractionation caused by gaseous diffusion
        b = 27.0  # effective fractionation caused by carboxylating enzyme
        ci_over_ca = (delta - a) / (b - a)

        return ci_over_ca

    def delta_c4(self, delta, temp):
        """
        Carbon isotope discrimination (delta) of leaf dry matter

        Parameters:
        -----------
        delta : float
            Carbon isotope discrimination
        temp : float
            growing-season mean temperature. mGDD0 is defined as the annual sum of
            temperatures above 0degC (GDD - growing degree days) divided by the
            length of the period with temperatures above 0degC

        Returns:
        --------
        ci_over_ca : float


        References:
        ----------
        * Henderson SA et al. (1992) Short-term measurements of carbon isotope
          discrimination in several C4 species. Functional Plant
          Biology, 19:263-285.
        * von Caemmerer S et al. (2014) Carbon isotope discrimination as a tool
          to explore C4 photosynthesis, Journal of Experimental Botany,
          doi:10.1093/jxb/eru127.

        """

        # the fractionation during diffusion of CO2 in air
        a = 4.4

        #  is the fractionation during Rubisco carboxylation
        b3 = 30.0

        # fractionation associated with PEP carboxylation
        # and the preceding isotopic equilibrium during dissolution
        # and conversion to bicarbonate
        b4 = -(9.483 * 1000.0) / (273.0 + temp) + 23.89 + 2.2

        # fractionation during the leakage of CO2 out of the bundle sheath
        # cells
        s = 1.8

        # leakiness
        phi = 0.21

        #delta = a + (b4 + phi * (b3 - s) - a) * ci_over_ca
        ci_over_ca = (delta - a) / (b4 + phi * (b3 - s) - a)

        return ci_over_ca

    def calculate_g1(self, ci_ca, vpd):
        return (ci_ca * np.sqrt(vpd)) / (1.0 - ci_ca)

if __name__ == "__main__":

    L = FitLeafIsotope(fdir="data/raw_data/isotope_data/",
                       ofdir="data/processed/",
                       fname="global_13C_dataset.csv",
                       pft_fname="global_13C_dataset_PFTs_from_wang_han_MDK_fixes.csv",
                       vpd_fname="CRU1.0_GDD-VPD.csv",
                       ofname="g1_isotope_screened.csv")

    L.main()
