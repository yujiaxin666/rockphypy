:py:mod:`rockphypy.Perm`
========================

.. py:module:: rockphypy.Perm

.. autoapi-nested-parse::

   Recreation and modifcations of the Permeabilities models in Rock physics handbook matlab tools.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.Perm.Permeability




.. py:class:: Permeability

   
   Different permeability models.
















   ..
       !! processed by numpydoc !!
   .. py:method:: Kozeny_Carman(phi, d)
      :staticmethod:

      
      Describe the permeability in a porous medium using Kozeny-Carman equation assuming the turtuosity tau=sqrt(2), 1/B=2.5 for unconsolidated monomodal sphere pack.

      :param phi: porosity
      :type phi: float or array-like
      :param d: pore diameter in m.
      :type d: float

      .. rubric:: Examples

      >>> phi= np.linspace(0.01,0.35,100)
      >>> d= 250
      >>> k= Kozeny_Carman(phi, d)
      >>> plt.semilogy(phi, k )

      :returns: *float or array-like* -- k (m^2): the resulting permeability is in the same units as d^2















      ..
          !! processed by numpydoc !!

   .. py:method:: Kozeny_Carman_Percolation(phi, phic, d, B)
      :staticmethod:

      
      The Kozeny−Carman relations incorporating the percolation effect

      :param phi: porosity
      :type phi: float or array-like
      :param phic: percolation porosity
      :type phic: float
      :param d: pore diameter
      :type d: float
      :param B: geometric factor that partly accounts for the irregularities of pore shapes.
      :type B: float

      :returns: *float or array-like* -- k (m^2): the resulting permeability is in the same units as d^2















      ..
          !! processed by numpydoc !!

   .. py:method:: Owolabi(phi, Swi)
      :staticmethod:

      
      Estimate the permeability in uncosonlidated sands of Pleistocene to Oligocene age in Eastern Niger Delta from log derived porosityand irreducible water saturation.

      :param phi: porosity
      :type phi: float or array-like
      :param Swi: irreducible water-saturation from welllogs
      :type Swi: float or array-like

      :returns: *float or array-like* -- k_oil, k_gas: permeabilities in mD for oil and gas sand reservoir, respectively















      ..
          !! processed by numpydoc !!

   .. py:method:: Perm_logs(phi, Swi)
      :staticmethod:

      
      Various empirical correlations of between permeability, porosity and irreducible water-saturation from welllogs. Models includs Tixier, Timur, Coates and Coates-Dumanoir.

      :param phi: porosity
      :type phi: float or array-like
      :param Swi: irreducible water-saturation from welllogs
      :type Swi: float or array-like

      :returns: * *float or array-like* -- k_tixier, k_Timur , k_coates, k_coates_Dumanoir: different permeability estimations, in the unit of mD
                * *Assumptions*
                * *-----------*
                * *- The functional forms used in these equations have to be calibrated, whenever possible, to site-specific data.*
                * *- The rock is isotropic.*
                * *- Fluid-bearing rock is completely saturated.*















      ..
          !! processed by numpydoc !!

   .. py:method:: Panda_Lake(d, C, S, tau, phi)
      :staticmethod:

      
      Modified Kozeny-carman relation incorpating the contribution of grain size variation and sorting using Manmath N. Panda and Larry W. Lake relation.

      :param d: mean particles size in um.
      :type d: float
      :param C: coefficient of variation of particles size distribution
      :type C: float
      :param S: skewness of particles size distribution
      :type S: float
      :param tau: tortuosity factor
      :type tau: float
      :param phi: porosity
      :type phi: float or array-like

      :returns: *float or array-like* -- k (md): permeability

      .. rubric:: References

      - Estimation of Single-Phase permeability from parameters of particle-Size Distribution, Manmath N. Panda and Larry W. Lake, AAPG 1994.















      ..
          !! processed by numpydoc !!

   .. py:method:: Panda_Lake_cem(phi, d)
      :staticmethod:

      
      Quantify the effects of cements on the single phase permeability estimate of unconsolidated sand using Panda & Lake model

      :param phi: porosity
      :type phi: float or array-like
      :param d: mean particles size in um
      :type d: float

      :returns: *float or array-like* -- k (md): permeability















      ..
          !! processed by numpydoc !!

   .. py:method:: Revil(phi, d)
      :staticmethod:

      
      Estimate permeability in very shaly rock using Revil et al. 1997

      :param phi: porosity
      :type phi: float or array-like
      :param d: mean particles size in um
      :type d: float

      :returns: *float or array-like* -- k (md): permeability















      ..
          !! processed by numpydoc !!

   .. py:method:: Fredrich(phi, d, b)
      :staticmethod:

      
      Compute permability considering Pore Geometry and Transport Properties of Fontainebleau Sandstone

      :param phi: porosity>10%
      :type phi: float or array-like
      :param d: _description_
      :type d: float
      :param b: shape factor b is equal to 2 for circular tubes and equal to 3 for cracks.
      :type b: float

      :returns: * *float or array-like* -- k (md): permeability
                * *References* -- ----------
                * *- Fredrich, J. T., Greaves, K. H., & Martin, J. W. (1993, December). Pore geometry and transport properties of Fontainebleau sandstone. In International journal of rock mechanics and mining sciences & geomechanics abstracts (Vol. 30, No. 7, pp. 691-697). Pergamon.*















      ..
          !! processed by numpydoc !!

   .. py:method:: Bloch(S, C, D)
      :staticmethod:

      
      Predict porosity and permeability in sandstones prior to drilling using Bloch empirical relations obtain in Yacheng field.

      :param S: Trask sorting coefficient
      :type S: float
      :param C: Rigid grain content in frac
      :type C: float
      :param D: Grain size in mm
      :type D: float

      :returns: *float or array-like* -- phi, k: porosity (frac) and permeability (mD), respectively















      ..
          !! processed by numpydoc !!

   .. py:method:: Bernabe(phi, crf, w, r)
      :staticmethod:

      
      Bernabe models permit to compute the permeability and porosity of strongly pressure dependent pores such as cracks and approximately constant pores associated with tubes and nodal pores.

      :param phi: total porosity
      :type phi: float or array-like
      :param crf: crack fraction in pore volume
      :type crf: float
      :param w: width or aperture of the equivalent crack in um
      :type w: float
      :param r: radius of the tube in um
      :type r: float

      :returns: * *float or array-like* -- k (md): total permeability
                * *References* -- ----------
                * *- Bernabe, Y. (1991). Pore geometry and pressure dependence of the transport properties in sandstones. Geophysics, 56(4), 436-446.*















      ..
          !! processed by numpydoc !!


