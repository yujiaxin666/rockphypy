:py:mod:`rockphypy.Emp`
=======================

.. py:module:: rockphypy.Emp


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.Emp.Empirical




.. py:class:: Empirical

   
   Empirical relations that widely applied
















   ..
       !! processed by numpydoc !!
   .. py:method:: krief(phi, Kg, Gg)
      :staticmethod:

      
      Compute porous background elastic constants as a function of porosity according to Krief model.

      :param phi: porosity in porous rock
      :type phi: float or array
      :param Kg: grain bulk modulus
      :type Kg: float
      :param Gg: grain shear modulus
      :type Gg: float

      :returns: *float or array* -- lamda,G,K















      ..
          !! processed by numpydoc !!

   .. py:method:: esti_VS(Vp, vsh)
      :staticmethod:

      
      Estimate, using the Greenberg-Castagna empirical relations, the shearwave velocity in a brine-saturated shaly sandstone with vp
      since we only assume two minearl phases: so L=2, X1= 1-vsh, X2= vsh

      :param Vp: compressional velocities, m/s
      :type Vp: float or array-like
      :param vsh: shale volume
      :type vsh: float or array-like

      :returns: * *float or array-like* -- vs (m/s): estimated shear wave velocities
                * *References* -- ----------
                * *- handbook of rock physics P516*















      ..
          !! processed by numpydoc !!

   .. py:method:: han(phi, C)
      :staticmethod:

      
      Han (1986) found empirical regressions relating ultrasonic (laboratory) velocities to porosity and clay content.effective pressure is 20Mpa

      :param phi: porosity
      :type phi: float or array-like
      :param C: clay volume fraction
      :type C: float or array-like

      :returns: *float or array-like* -- P and S wave velocities















      ..
          !! processed by numpydoc !!

   .. py:method:: ehrenberg(Z)
      :staticmethod:

      
      porosity reference trend for Norwegian Sea sandstone. Note that the functional form of the porosity model is not published in Ehrenberg (1990). It is obtained by linear regression of the digitized data point from the original plot in the paper.

      :param Z: burial depth below see floor in Km
      :type Z: float or array

      .. rubric:: References

      Ehrenberg, S., 1990, Relationship between diagenesis and reservoir quality in sandstones of the Garn Formation, Haltenbanken, mid-Norwegian continental shelf: AAPG bulletin, 74, no. 10, 1538

      :returns: *float or array* -- porosity















      ..
          !! processed by numpydoc !!

   .. py:method:: yu_segment_trend(Z)
      :staticmethod:

      
      Reference trend for Norwegian sea normally buried clean sandstones

      :param Z: burial depth below see floor in m
      :type Z: float or array

      :returns: *float or array* -- P wave velocities















      ..
          !! processed by numpydoc !!

   .. py:method:: ramm_porosity(Z, HB=True)
      :staticmethod:

      
      porosity reference trend according to Ramm & Bjørlykke (1994)

      :param Z: burial depth wrt. sea floor in m
      :type Z: float or array
      :param HB: if True: only show the regression line for halten bakken area porosity data False: The regression line for all porosity from north sea and norwegian sea, by default True
      :type HB: bool, optional

      :returns: *float or array* -- porosity















      ..
          !! processed by numpydoc !!

   .. py:method:: ramm_porosity_segment(Z)
      :staticmethod:

      
      segment porosity reference trend according to Ramm & Bjørlykke (1994) considering the mechanical and chemical compaction

      :param Z: burial depth wrt. sea floor in m
      :type Z: float or array

      :returns: *float or array* -- porosity















      ..
          !! processed by numpydoc !!

   .. py:method:: empirical_StPeter(Pe, sample=1)
      :staticmethod:

      
      compute the Vp and Vs for st peter sandstone using the empirical relationship in the form of V= A+KPe-Be^(-DPe)

      :param Pe: effective pressure in Kbar, 1kbar= 100Mpa
      :type Pe: float or array
      :param sample: 1-sample 1, phi= 0.205. 2-sample 2, phi= 0.187, by default 1
      :type sample: int, optional

      :returns: *float or array* -- Vp,Vs in km















      ..
          !! processed by numpydoc !!

   .. py:method:: Scherbaum(Z)
      :staticmethod:

      
      velocity depth trend for Lower and Middle Buntsandstein

      :param Z: burial depth wrt. sea floor in m
      :type Z: float or array

      .. rubric:: References

      Scherbaum, F., 1982. Seismic velocities in sedimentary rocks—indicators of subsidence and uplift?. Geologische Rundschau, 71(2), pp.519-536.

      :returns: *float or array* -- P wave velocities in m/s















      ..
          !! processed by numpydoc !!

   .. py:method:: Sclater(phi)
      :staticmethod:

      
      Sclater-Christie exponential curve for sandstone

      :param phi: porosity
      :type phi: float or array

      :returns: **Z** (*float or array*) -- depth wrt. sea floor in km.















      ..
          !! processed by numpydoc !!

   .. py:method:: Storvoll(Z)
      :staticmethod:

      
      Storvoll velocity compaction trend. The trend is for shale and shaly sediments but also used for siliciclastic rock like sandstone

      :param Z: depth wrt. sea floor in m.
      :type Z: float or array

      :returns: *float or array* -- Vp meters per second















      ..
          !! processed by numpydoc !!

   .. py:method:: Hillis(Z)
      :staticmethod:

      
      compaction trend for Bunter Sandstone in North sea

      :param Z: depth below sea bed (in kilometers)
      :type Z: float or array

      :returns: *float or array* -- Vp km/s















      ..
          !! processed by numpydoc !!

   .. py:method:: Japsen(Z)
      :staticmethod:

      
      a segmented linear velocity–depth function, These equations are considered as approximation for bunter sandstone trend although they are originally for bunter shale. proposed by Japsen 1999

      :param Z: depth below sea bed in m.
      :type Z: float or array

      :returns: *float or array* -- Vp m/s















      ..
          !! processed by numpydoc !!

   .. py:method:: hjelstuen(Z)
      :staticmethod:

      
      Velocity-depth relationships for the Bjørna-Sørkapp margin deposits. note:  the seismic velocities are not directly comparable with velocities from sonic logs (because of the different frequencies), and the velocity-depth profile of Hjelstuen et al. (1996) has not been corrected for uplift and erosion

      :param Z: Z< 3.8km
      :type Z: float or array

      :returns: *float or array* -- V: m/s















      ..
          !! processed by numpydoc !!

   .. py:method:: Cp(phi)
      :staticmethod:

      
      The coordination number n depends on porosity, as shown by
      Murphy 1982.

      :param phi: total porosity , for a porosity of 0.4, n=8.6
      :type phi: float or array

      :returns: *float or array* -- coordination number















      ..
          !! processed by numpydoc !!


