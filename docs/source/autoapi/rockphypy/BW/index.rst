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
   .. py:method:: dz_dp(T_pr)

      
      Values for dZ/dPpr obtained from equation 10b in Batzle and Wang (1992).
















      ..
          !! processed by numpydoc !!

   .. py:method:: pseudo_p_t(T, G)

      
      Calculate the pseudoreduced temperature and pressure according to Thomas et al. 1970.

      :param P: Pressure in MPa
      :type P: float or array-like
      :param T: Temperature in °C
      :type T: float or array-like
      :param G: Gas gravity
      :type G: float or array-like

      :returns: *float or array-like* -- Ta: absolute temperature
                Ppr:pseudoreduced pressure
                Tpr:pseudoreduced temperature















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_co2(T, G)

      
      Compute CO2 properties as a function of temperature and pressure using modified Batzle-Wang equations

      :param P: Pressure in MPa
      :type P: float or array-like
      :param T: Temperature in °C
      :type T: float or array-like
      :param G: Gas gravity
      :type G: float or array-like

      :returns: *float or array-like* -- rho (g/cc): gas density
                K (GPa): bulk modulus

      .. rubric:: References

      Xu, H. (2006). Calculation of CO2 acoustic properties using Batzle-Wang equations. Geophysics, 71(2), F21-F23.















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_gas(T, G)

      
      Estimate the Gas density and bulk modulus at specific temperature and pressure.

      :param P: Pressure in MPa
      :type P: float or array-like
      :param T: Temperature in °C
      :type T: float or array-like
      :param G: Gas gravity
      :type G: float or array-like

      :returns: *float or array-like* -- rho: Gas density
                K: Gas bulk modulus















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_oil(T, den)

      
      Estimate the oil density and bulk modulus at specific temperature and pressure.

      :param P: Pressure in MPa
      :type P: float or array-like
      :param T: Temperature in °C
      :type T: float or array-like
      :param den: oil density in g/cm3
      :type den: float or array-like

      :returns: *float or array-like* -- rho: oil density
                K: oil bulk modulus















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_go(T, den, G, Rg)

      
      compute density and bulk modulus of live oil.

      :param P: Pressure in MPa
      :type P: float or array-like
      :param T: Temperature in °C
      :type T: float or array-like
      :param den: oil density in g/cm3
      :type den: float or array-like
      :param G: gas gravity
      :type G: float or array-like
      :param Rg: the volume ratio of liberated gas to remaining oil at atmospheric pressure and 15.6°C, Liter/Liter
      :type Rg: float or array-like

      :returns: *float or array-like* -- v (m/s): velocity
                rho_g (g/cm3): true density of live oil at saturation
                K (GPa): true bulk modulus of live oil at saturation















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_water(P)

      
      Compute the density and bulk modulus of pure water as a function of temperature and pressure using Batzle and Wang (1992).

      :param T: Temperature in °C
      :type T: float or array-like
      :param P: Pressure in MPa
      :type P: float or array-like

      :returns: *float or array-like* -- rho_w (g/cm3): density of pure water















      ..
          !! processed by numpydoc !!

   .. py:method:: v_water(P)

      
      Acoustic velocity of pure water as a function of temperature
      and pressure using Batzle and Wang (1992).

      :param T: Temperature in °C
      :type T: float or array-like
      :param P: Pressure in MPa
      :type P: float or array-like

      :returns: *float or array-like* -- v_w (m/s): acoustic velocity of pure water















      ..
          !! processed by numpydoc !!

   .. py:method:: rho_K_brine(P, S)

      
      Calculation of the density and bulk modulus of brine (NaCl) as a function of temperature, salinity and pressure using Batzle and Wang (1992).

      :param T: Temperature in °C
      :type T: float or array-like
      :param P: Pressure in MPa
      :type P: float or array-like
      :param S: weight fraction of sodium chloride in ppm/1e6
      :type S: float or array-like

      :returns: *float or array-like* -- rho_b (g/cm3): the density of brine
                K_b (GPa):bulk modulus of brine















      ..
          !! processed by numpydoc !!

   .. py:method:: v_brine(P, S)

      
      Calculte the acoustic velocity of brine as a function of temperature, salinity and pressure using Batzle and Wang (1992).

      :param T: Temperature in °C
      :type T: float or array-like
      :param P: Pressure in MPa
      :type P: float or array-like
      :param S: weight fraction of sodium chloride in ppm/1e6
      :type S: float or array-like

      :returns: *float or array-like* -- v_b (m/s): the velocity of brine















      ..
          !! processed by numpydoc !!


