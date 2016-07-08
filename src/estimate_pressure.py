# -*- coding: UTF-8 -*-
"""
The hypsometric equation, also known as the thickness equation, relates an
atmospheric pressure ratio to the equivalent thickness of an atmospheric
layer under the assumptions of constant temperature and gravity.
It is derived from the hydrostatic equation and the ideal gas law.
"""
from numpy import log, exp
import numpy as np
import sys

__author__  = "Martin De Kauwe"
__version__ = "1.0 (18.02.2016)"
__email__   = "mdekauwe@gmail.com"


def estimate_pressure(tair, elev):
    """
    Estimate pressure using the hypsometric equation, An equation relating
    the thickness, h, between two isobaric surfaces to the mean virtual
    temperature of the layer.

    Parameters:
    -----------
    tair : float
        air temperature (deg C)
    elev : float
        elevation (m)

    Returns:
    -----------
    pressure : float
        air pressure [Pa]

    References:
    ----------
    * http://glossary.ametsoc.org/wiki/Hypsometric_equation
    * Jurgen sent me this function
    """
    C_TO_K = 273.15
    RGAS = 8.314                # J mol-1 K-1
    MILLIBAR_2_PA = 100.

    mv = 18.0     # molecular weight of water vapor (g per mol)
    ma = 28.9     # molecular weight of dry air (g per mol)
    g = 9.8       # m s-2
    qmol = 0.622

    Tk = tair + C_TO_K
    # first approx the pressure with the air temperature (mb)
    press = 1013.0 / (exp(g / (RGAS / (ma / 1000.0)) / Tk * elev))
    esat = 6.112 * exp(17.61 * (Tk - 273.15) / (Tk - 29.65)) # mb

    # g water vapor kg-1 moist air
    qsat = esat * 0.622 / (press - 0.378 * esat) / 1000.0
    qs = qmol * mv / ma  # g/kg
    ea = qs / 0.622  # in mb

    qs_kg = qs / 1000.0 # kg H2O vapor / kg air
    ea_pa = ea * 100.0  # in Pa

    # Estimate air pressure from hypsometric adjustment for elevation (mb)
    TvK = Tk * (1.0 + 0.61 * qs_kg)
    press = 1013.0 / (exp(g / (RGAS / (ma / 1000.0)) / TvK * elev))
    press *= MILLIBAR_2_PA

    return (press)
