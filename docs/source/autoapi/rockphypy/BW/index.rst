:py:mod:`rockphypy.BW`
======================

.. py:module:: rockphypy.BW

.. autoapi-nested-parse::

   Batzle and Wang functionalities

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.BW.BW




.. py:class:: BW

   
   Effective CO2, natural gas, brine and oil property calculation using original and modified Batzle-Wang equations.
















   ..
       !! processed by numpydoc !!
   .. py:method:: dz_dp(P_pr, T_pr)
      :staticmethod:

      
      Values for dZ/dPpr obtained from equation 10b in Batzle and Wang (1992).
















      ..
          !! processed by numpydoc !!

   .. py:method:: pseudo_p_t(P, T, G)
      :staticmethod:

      
      Calculate the pseudoreduced temperature and pressure according to Thomas et al. 1970.

      :param P: Pressure in MPa
      :type P: float or array-like
      :param T: Temperature in °C
      :type T: float or array-like
      :param G: Gas gravity
      :type G: float

      :returns: *float or array-like* -- Ta: absolute temperature
                Ppr:pseudoreduced pressure
                Tpr:pseudoreduced temperature















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_co2(P, T, G)
      :staticmethod:

      
      Compute CO2 properties as a function of temperature and pressure using modified Batzle-Wang equations

      :param P: Pressure in MPa
      :type P: float or array-like
      :param T: Temperature in °C
      :type T: float or array-like
      :param G: Gas gravity
      :type G: float

      :returns: *float or array-like* -- rho (g/cc): gas density
                K (GPa): bulk modulus

      .. rubric:: References

      Xu, H. (2006). Calculation of CO2 acoustic properties using Batzle-Wang equations. Geophysics, 71(2), F21-F23.















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_gas(P, T, G)
      :staticmethod:

      
      Estimate the Gas density and bulk modulus at specific temperature and pressure.

      :param P: Pressure in MPa
      :type P: float or array-like
      :param T: Temperature in °C
      :type T: float or array-like
      :param G: Gas gravity
      :type G: float

      :returns: *float or array-like* -- rho: Gas density (g/cm3)
                K: Gas bulk modulus (GPa)















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_oil(P, T, den)
      :staticmethod:

      
      Estimate the oil density and bulk modulus at specific temperature and pressure.

      :param P: Pressure in MPa
      :type P: float or array-like
      :param T: Temperature in °C
      :type T: float or array-like
      :param den: oil density in g/cm3
      :type den: float

      :returns: *float or array-like* -- rho: oil density (g/cm3)
                K: oil bulk modulus (GPa)















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_go(P, T, den, G, Rg)
      :staticmethod:

      
      compute density and bulk modulus of live oil.

      :param P: Pressure in MPa
      :type P: float or array-like
      :param T: Temperature in °C
      :type T: float or array-like
      :param den: oil density in g/cm3
      :type den: float
      :param G: gas gravity
      :type G: float
      :param Rg: the volume ratio of liberated gas to remaining oil at atmospheric pressure and 15.6°C, Liter/Liter
      :type Rg: float

      :returns: *float or array-like* -- rho_g (g/cm3): true density of live oil at saturation
                K (GPa): true bulk modulus of live oil at saturation















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_water(T, P)
      :staticmethod:

      
      Compute the density and bulk modulus of pure water as a function of temperature and pressure using Batzle and Wang (1992).

      :param T: Temperature in °C
      :type T: float or array-like
      :param P: Pressure in MPa
      :type P: float or array-like

      :returns: *float or array-like* -- rho_w (g/cm3): density of pure water
                K_w (Gpa): bulk modulus of pure water















      ..
          !! processed by numpydoc !!

   .. py:method:: v_water(T, P)
      :staticmethod:

      
      Acoustic velocity of pure water as a function of temperature
      and pressure using Batzle and Wang (1992).

      :param T: Temperature in °C
      :type T: float or array-like
      :param P: Pressure in MPa
      :type P: float or array-like

      :returns: *float or array-like* -- v_w (m/s): acoustic velocity of pure water















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_brine(T, P, S)
      :staticmethod:

      
      Calculation of the density and bulk modulus of brine (NaCl) as a function of temperature, salinity and pressure using Batzle and Wang (1992).

      :param T: Temperature in °C
      :type T: float or array-like
      :param P: Pressure in MPa
      :type P: float or array-like
      :param S: weight fraction of sodium chloride in ppm/1e6
      :type S: float

      :returns: *float or array-like* -- rho_b (g/cm3): the density of brine
                K_b (GPa):bulk modulus of brine















      ..
          !! processed by numpydoc !!

   .. py:method:: v_brine(T, P, S)
      :staticmethod:

      
      Calculte the acoustic velocity of brine as a function of temperature, salinity and pressure using Batzle and Wang (1992).

      :param T: Temperature in °C
      :type T: float or array-like
      :param P: Pressure in MPa
      :type P: float or array-like
      :param S: weight fraction of sodium chloride in ppm/1e6
      :type S: float

      :returns: *float or array-like* -- v_b (m/s): the velocity of brine















      ..
          !! processed by numpydoc !!

   .. py:method:: co2_brine(temperature, pressure, salinity, Sco2, brie_component=None, bw=False)
      :staticmethod:

      
      compute the effective properties of critical Co2 brine mixture depending on temperature, pressure and salinity of the brine, as well as the saturation state.

      :param temperature:
      :type temperature: degree
      :param pressure: pore pressure, not effective stress
      :type pressure: Mpa
      :param salinity: The weight fraction of NaCl, e.g. 35e-3
                       for 35 parts per thousand, or 3.5% (the salinity of
                       seawater).
      :type salinity: ppm
      :param Sco2: Co2 saturation
      :type Sco2: frac
      :param brie_component: if None: uniform saturation. otherwise patchy saturation according to brie mixing
      :type brie_component: num

      :returns: *den_mix (g/cc)* -- mixture density
                Kf_mix (GPa): bulk modulus















      ..
          !! processed by numpydoc !!


