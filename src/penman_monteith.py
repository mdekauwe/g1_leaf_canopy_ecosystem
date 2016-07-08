# -*- coding: UTF-8 -*-

from numpy import log, exp
import numpy as np
import sys

__author__  = "Martin De Kauwe"
__version__ = "1.0 (18.02.2016)"
__email__   = "mdekauwe@gmail.com"

class PenmanMonteith(object):

    """
    Penman-Monteith equation to calculate canopy transpiration.

    - Class also contain a method to invert canopy conducance (gc) if
      transpiration is already known

    References:
    -----------
    * Monteith and Unsworth (1990) Principles of Environmental
      Physics, pg. 247. Although I have removed the soil heat flux as G'DAY
      calculates soil evaporation seperately.

    """

    def __init__(self, dz0v_dh=0.075, z0h_z0m=0.1, d=0.67, use_ustar=False):

        """
        Parameters:
        -----------
        cp : float
            specific heat of dry air [MJ kg-1 degC-1]
        vk : float
            von Karman's constant [unitless]
        epsilon : float
            ratio molecular weight of water vap/dry air
        zele_sea : float
            elevation above sea level [m]
        dz0v_dh : float
            rate change of roughness for momentum with height
        displace_ratio : float
            zero plain displacement height
        z0h_z0m : float
            Ratio of the roughness length for heat to the roughness length for
            momentum, see comment in method below!!!
        """

        self.CP = 1010.0                 # specific heat of dry air (j kg-1 k-1)
        self.VK = 0.41                   # von Karman constan
        self.J_TO_MJ = 1.0E-6
        self.C_TO_K = 273.15
        self.dz0v_dh = dz0v_dh
        self.displace_ratio = d          # zero plan displacement height
        self.z0h_z0m = z0h_z0m
        self.RGAS = 8.314                # J mol-1 K-1
        self.H2OLV0 = 2.501E6            # latent heat H2O (J kg-1)
        self.H2OMW = 18E-3               # mol mass H20 (kg mol-1)
        self.MASS_AIR = 29.0E-3          # mol mass air (kg mol-1)
        self.use_ustar = use_ustar       # calc ga using ustar

    def calc_evaporation(self, vpd, wind, rnet, tair, press, gs, canht=None,
                         ustar=None, G=None):
        """
        Parameters:
        -----------
        vpd : float
            vapour pressure def [Pa]
        wind : float
            average daytime wind speed [m s-1]
        rnet : float
            net radiation [W m-2]
        tair : float
            temperature [degC]
        press : float
            average daytime pressure [Pa]
        gs : float
            stomatal conductance [mol m-2 s-1]
        canht : float
            canopy height [m]
        ustar : float
            friction velocity [m s-1]
        G : float
            soil heat flux [W m-2]

        Returns:
        --------
        trans : float
            transpiration [mol H20 m-2 s-1]
        """

        # use friction velocity
        if self.use_ustar:
            ga = self.calc_bdary_layer_conduc_from_ustar(wind, ustar, press,
                                                         tair)
        else:
            ga = self.canopy_bdary_layer_conduct(canht, wind, press, tair)
        lambdax = self.calc_latent_heat_of_vapourisation(tair)
        gamma = self.calc_pyschrometric_constant(lambdax, press)
        slope = self.calc_slope_of_sat_vapour_pressure_curve(tair)

        # Total leaf conductance to water vapour
        #gv = 1.0 / (1.0 / gs + 1.0 / ga)

        if G is None:
            G = 0.0 # ground heat flux

        arg1 = slope * (rnet - G) + ga * self.MASS_AIR * self.CP * vpd
        arg2 = slope + gamma * (1.0 + ga / gs)
        LE = arg1 / arg2 # W m-2

        # mol H20 m-2 s-1
        transpiration = LE / lambdax
        transpiration = np.where(transpiration < 0.0, 0.0, transpiration)

        return transpiration

    def invert_penman(self, vpd, wind, rnet, tair, press, trans, canht=None,
                      ustar=None, G=None):
        """
        Invert Penman-Monteith eqn to obtain canopy conductance

        Parameters:
        -----------
        vpd : float
            vapour pressure def [Pa]
        wind : float
            average daytime wind speed [m s-1]
        rnet : float
            net radiation [W m-2]
        tair : float
            temperature [degC]
        press : float
            average daytime pressure [Pa]
        trans : float
            transpiration [mol H20 m-2 s-1]
        canht : float
            canopy height [m]
        ustar : float
            friction velocity [m s-1]
        G : float
            soil heat flux [W m-2]

        Returns:
        --------
        gc : float
            canopy conductance [mol m-2 s-1]

        Reference:
        ---------
        * Landsberg and Sands, eqn 2.53
        """
        if self.use_ustar:
            ga = self.calc_bdary_layer_conduc_from_ustar(wind, ustar, press,
                                                         tair)
        else:
            ga = self.canopy_bdary_layer_conduct(canht, wind, press, tair)

        lambdax = self.calc_latent_heat_of_vapourisation(tair)
        gamma = self.calc_pyschrometric_constant(lambdax, press)
        slope = self.calc_slope_of_sat_vapour_pressure_curve(tair)
        lambda_E = trans * lambdax

        if G is None:
            G = 0.0 # ground heat flux

        arg1 = ga * gamma * lambda_E
        arg2 = slope * (rnet - G) - (slope + gamma) * lambda_E
        arg3 = ga * self.MASS_AIR * self.CP * vpd

        return arg1 / (arg2 + arg3)

    def calc_decoupling_coefficent(self, wind, tair, press, gs, canht=None,
                                  ustar=None):
        """
        Calculate decoupling coefficient.
            - As omega -> 0, leaf surface becomes strongly coupled to the
              atmosphere.
            - As omega -> 1, leaf surfaces are poorly coupled to the
              atmosphere.

        Parameters:
        -----------
        wind : float
            average daytime wind speed [m s-1]
        tair : float
            temperature [degC]
        press : float
            average daytime pressure [Pa]
        gs : float
            stomatal conductance [mol H20 m-2 s-1]
        canht : float
            canopy height [m]
        ustar : float
            friction velocity [m s-1]

        Returns:
        --------
        omega : float
            decoupling coefficient (-)


        References:
        -----------
        * McNaughton and Jarvis 1986
        """
        if self.use_ustar:
            ga = self.calc_bdary_layer_conduc_from_ustar(wind, ustar, press,
                                                         tair)
        else:
            ga = self.canopy_bdary_layer_conduct(canht, wind, press, tair)

        lambdax = self.calc_latent_heat_of_vapourisation(tair)
        gamma = self.calc_pyschrometric_constant(lambdax, press)
        slope = self.calc_slope_of_sat_vapour_pressure_curve(tair)

        epsilon = slope / gamma
        omega = (1.0 + epsilon) / (1.0 + epsilon + ga / gs)

        return (omega)

    def calc_bdary_layer_conduc_from_ustar(self, wind, ustar, press, tair):
        """
        Calculate boundary layer conductance using measured friction velocity,
        ustar

        Parameters:
        -----------
        wind : float
            average daytime wind speed [m s-1]
        ustar : float
            friction velocity [m s-1]

        Returns:
        --------
        ga : float
            canopy boundary layer conductance (mol m-2 s-1)

        References:
        ----------
        * Monteith & Unsworth p. 341 eqn. 17.8
        """

        # Convert from m s-1 to mol m-2 s-1
        # - note conversion in Jones '92 is mmol to mmol, but units cancel
        Tk = tair + self.C_TO_K
        cmolar = press / (self.RGAS * Tk)

        ga = 1.0 / (wind / ustar**2 + 6.2 * ustar**-0.667)
        ga *= cmolar

        return (ga)

    def canopy_bdary_layer_conduct(self, canht, wind, press, tair):
        """  Canopy boundary layer conductance, ga (from Jones 1992 p 68)

        Notes:
        ------
        'Estimates of ga for pine canopies from LAI of 3 to 6 vary from
        3.5 to 1.1 mol m-2 s-1  (Kelliher et al., 1993; Juang et al., 2007).'
        Drake et al, 2010, 17, pg. 1526.

        References:
        ------------
        * Jones 1992, pg. 67-8.
        * Monteith and Unsworth (1990), pg. 248. Note this in the inverted form
          of what is in Monteith (ga = 1 / ra)
        * Allen et al. (1989) pg. 651.
        * Gash et al. (1999) Ag forest met, 94, 149-158.

        Parameters:
        -----------
        wind : float
            average daytime wind speed [m s-1]
        press : float
            atmospheric pressure (Pa)
        tair : float
            air temperature (deg C)
        canht : float
            canopy height (m)

        Returns:
        --------
        ga : float
            canopy boundary layer conductance [mol m-2 s-1]
        """

        # Convert from m s-1 to mol m-2 s-1
        # - note conversion in Jones '92 is mmol to mmol, but units cancel
        Tk = tair + self.C_TO_K
        cmolar = press / (self.RGAS * Tk)

        # roughness length for momentum
        z0m = self.dz0v_dh * canht

        # z0h roughness length governing transfer of heat and vapour [m]
        # *Heat tranfer typically less efficent than momentum transfer. There is
        #  a lot of variability in values quoted for the ratio of these two...
        #  JULES uses 0.1, Campbell and Norman '98 say z0h = z0m / 5. Garratt
        #  and Hicks, 1973/ Stewart et al '94 say z0h = z0m / 7. Therefore for
        #  the default I am following Monteith and Unsworth, by setting the
        #  ratio to be 1, the code below is identical to that on page 249,
        #  eqn 15.7
        z0h = self.z0h_z0m * z0m

        # zero plan displacement height [m]
        d = self.displace_ratio * canht

        arg1 = self.VK**2 * wind
        arg2 = log((canht - d) / z0m)
        arg3 = log((canht - d) / z0h)

        ga = (arg1 / (arg2 * arg3)) * cmolar

        return (ga)

    def calc_slope_of_sat_vapour_pressure_curve(self, tair):
        """ Constant slope in Penman-Monteith equation

        Parameters:
        -----------
        tavg : float
            average daytime temperature

        Returns:
        --------
        slope : float
            slope of saturation vapour pressure curve [Pa K-1]
        """
        arg1 = self.calc_sat_water_vapour_press(tair+0.1)
        arg2 = self.calc_sat_water_vapour_press(tair)
        slope = (arg1 - arg2) / 0.1

        return (slope)

    def calc_sat_water_vapour_press(self, tac):
        """ Calculate saturated water vapour pressure at temperature TAC

        Parameters:
        -----------
        tac : float
            Celsius

        Returns:
        --------
        sat : float
            units: Pa

        References:
        -----------
        * Jones 1992 p 110 (note error in a - wrong units)
        """

        return (613.75 * exp(17.502 * tac / (240.97 + tac)));

    def calc_pyschrometric_constant(self, lambdax, press):
        """ Psychrometric constant ratio of specific heat of moist air at
        a constant pressure to latent heat of vaporisation.

        Parameters:
        -----------
        press : float
            air pressure (Pa)
        lambda : float
             latent heat of water vaporization (J mol-1)

        Returns:
        --------
        gamma : float
            pyschrometric constant [Pa K-1]
        """

        return (self.CP * self.MASS_AIR * press / lambdax)

    def calc_latent_heat_of_vapourisation(self, tair):
        """
        Latent heat of water vapour at air temperature

        Returns:
        -----------
        lambda : float
            latent heat of water vaporization [J mol-1]
        """

        return ((self.H2OLV0 - 2.365E3 * tair) * self.H2OMW)



if __name__ == "__main__":

    KPA_2_PA = 1000.0

    # these numbers are entirely made up!
    tair = 25.0
    vpd = 1.5 * KPA_2_PA
    wind = 1.5
    rnet = 12.0 # W m-2
    press = 101.325 * KPA_2_PA
    gs = 0.14915878
    canht = 12.0

    # Coniferous forest, based on Jones 1976, from Jones '92, pg 67
    P = PenmanMonteith(dz0v_dh=0.075, z0h_z0m=0.1)
    trans = P.calc_evaporation(vpd, wind, rnet, tair, press, gs, canht)
    # mmol m-2 s-1
    print trans * 1000.0
    print gs

    gc = P.invert_penman(vpd, wind, rnet, tair, press, trans, canht)
    print gc
    print
