#!/usr/bin/env python

"""
Fit Belinda's stomatal conductance model

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


class FitMedlyn(object):
    """
    Fit Belinda's stomatal conductance model using the lmfit package

    Reference:
    ---------
    * Medlyn et al. (2011) Global Change Biology (2011) 17, 2134-2144.

    """

    def __init__(self, fluxnet=True):
        self.report_fname = None
        self.param_fits = ["g0", "g1"]
        self.model_fit = None
        self.fluxnet = fluxnet

    def gs_model(self, vpd, gpp, co2, g0, g1):
        return g0 + 1.6 * (1.0 + (g1 / np.sqrt(vpd))) * (gpp / co2)

    def residual(self, params, df, obs):
        g0 = params['g0'].value
        g1 = params['g1'].value
        if self.fluxnet:
            model = self.gs_model(df["VPD_f"], df["GPP_f"], df["CO2"], g0, g1)
        else:
            model = self.gs_model(df["VPD"], df["Photo"], df["CO2S"], g0, g1)
        return (obs - model)

    def setup_model_params(self):
        """ Setup parameters """
        params = Parameters()
        params.add('g0', value=0.0, vary=False)
        params.add('g1', value=2.0, min=0.0)

        return params

    def print_fit_to_screen(self, result):
        for name, par in result.params.items():
            print '%s = %.4f +/- %.4f ' % (name, par.value, par.stderr)
        print

    def minimise_params(self, params, data, obs):
        try:
            result = minimize(self.residual, params, args=(data, obs))
            success = True
        except:
            result = -999
            success = False

        if (result.success == False or
            result.ier > 4 or
            result.errorbars == False):
            success = False

        return (result, success)



if __name__ == "__main__":

    main()
