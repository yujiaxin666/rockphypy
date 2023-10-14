:py:mod:`rockphypy`
===================

.. py:module:: rockphypy


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   AVO/index.rst
   Anisotropy/index.rst
   BW/index.rst
   EM/index.rst
   Emp/index.rst
   Fluid/index.rst
   GM/index.rst
   Perm/index.rst
   QI/index.rst
   utils/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.Anisotropy
   rockphypy.AVO
   rockphypy.BW
   rockphypy.Permeability
   rockphypy.GM
   rockphypy.Fluid
   rockphypy.EM
   rockphypy.Empirical
   rockphypy.QI
   rockphypy.utils




.. py:class:: Anisotropy

   
   Effective models, coordinate transform, anisotropic parameters and phase velocities that can be applied to anisotropic media.
















   ..
       !! processed by numpydoc !!
   .. py:method:: Thomsen(C11, C33, C13, C44, C66, den, theta)
      :staticmethod:

      
      Compute thomsen parameters and three phase velocities for weak anisotropic TI media with vertical symmetry axis.

      :param Cij: Stiffnesses in GPa
      :type Cij: float
      :param den: density of the effective medium, for backus average, the effective density can be computed using VRH method for Voigt average.
      :type den: float
      :param theta: angle of incidence
      :type theta: float or array-like

      :returns: *float or array-like* -- VP, VSV, VSH: wave velocities propagating along given direction















      ..
          !! processed by numpydoc !!

   .. py:method:: Thomsen_Tsvankin(C11, C22, C33, C12, C13, C23, C44, C55, C66)
      :staticmethod:

      
      Elastic constants of an orthorhombic elastic medium defined by Tsvankin’s notation for weak elastic anisotropy assuming the vertical symmetry axis is along the x3 direction.

      :param Cij: Stiffnesses in GPa
      :type Cij: float

      :returns: *floats* -- Thomsen-Tsvankin parameters















      ..
          !! processed by numpydoc !!

   .. py:method:: Backus(V, lamda, G)
      :staticmethod:

      
      Computes stiffnesses of a layered medium using backus average model.

      :param V: volumetric fractions of N isotropic layering materials
      :type V: float or array-like
      :param lamda: Lamé coefficients of N isotropic layering materials
      :type lamda: float or array-like
      :param G: shear moduli of N isotropic layering materials
      :type G: float or array-like

      :returns: *float or array-like* -- C11,C33,C13,C44,C66:Elastic moduli of the anisotropic layered media















      ..
          !! processed by numpydoc !!

   .. py:method:: Backus_log(Vp, Vs, Den, Depth)
      :staticmethod:

      
      Computes Backus Average from log data, notice that the Depth is 1d Vector including each top depth of layer and also the bottom of last layer.

      :param Vp: P wave velocities of layers [Vp1,Vp2...Vpn], Km/s, size N
      :type Vp: array
      :param Vs: S wave velocities of layers [Vs1,Vs2...Vsn],Km/s size N
      :type Vs: array
      :param Den: Densities of layers, size N
      :type Den: array
      :param Depth: 1d depth, ATTENTION: each depth point corresponds to the top of thin isotropic layer, the bottom of the sedimentary package is the last depth point. [dep1,dep2,,,,,,depn, depn+1], size N+1
      :type Depth: array

      :returns: *array-like* -- Stiffness coeffs and averaged density















      ..
          !! processed by numpydoc !!

   .. py:method:: vel_azi_HTI(C, Den, azimuth)
      :staticmethod:

      
      Given stiffnesses and density of the HTI medium, compute the azimuth dependent phase velocities.

      :param C: stiffness matrix of the HTI medium
      :type C: 2d array
      :param Den: density of the fractured medium
      :type Den: float
      :param azimuth: azimuth angle, degree
      :type azimuth: float or array like

      :returns: *float or array like* -- VP,VSH, VSV: phase velocities















      ..
          !! processed by numpydoc !!

   .. py:method:: vel_azi_VTI(C, Den, azimuth)
      :staticmethod:

      
      Given stiffnesses and density of the VTI medium, compute the azimuth dependent phase velocities.

      :param C: stiffness matrix of the VTI medium
      :type C: 2d array
      :param Den: density of the fractured medium
      :type Den: float
      :param azimuth: azimuth angle, degree
      :type azimuth: float or array like

      :returns: *float or array like* -- VP,VSH, VSV: phase velocities















      ..
          !! processed by numpydoc !!

   .. py:method:: Bond_trans(C, theta, axis=3)
      :staticmethod:

      
      Coordinate Transformations for stiffness matrix in 6x6 Voigt notation using Bond transformation matrix.

      :param C: original stiffness matrix
      :type C: 2d array
      :param theta: rotational angle
      :type theta: float
      :param axis: axis=1: fix 1-axis, rotate 2 and 3 axis, examples can be a TTI(Tilted TI) resulted from the rotation of VTI with horizontal aligned fracture sets wrt the vertical x3 axis. In this case, the input C should be a VTI matrix
                   axis=3: fix 3-axis, rotate 1 and 2 axis, E.g. seismic measurements of HTI media e.g caused by vertically aligned fractures. The angle theta may be assigned to be the angle between the fracture normal and a seismic line.
      :type axis: int, optional

      :returns: *2d array* -- C_trans, new stiffness matrix wrt to the original right-hand rectangular Cartesian coordinates

      .. rubric:: References

      - Bond, W., Jan. 1943, The mathematics of the physical properties of crystals, The Bell System Technical Journal, 1-72.















      ..
          !! processed by numpydoc !!


.. py:class:: AVO

   
   Exact and approximations of reflectivity in isotropic and anisotropic media.
















   ..
       !! processed by numpydoc !!
   .. py:method:: AVO_HTI(D1, D2, C1, C2, theta, azimuth)
      :staticmethod:

      
      Compute azimuth dependent PP reflectivity for wealy anisotropic HTI media using Ruger's approximation

      :param D1: density of the upper medium g/cm3
      :type D1: float or array-like
      :param D2: density of the lower medium g/cm3
      :type D2: float or array-like
      :param C1: stiffness matrix of the upper medium
      :type C1: 2D array
      :param C2: stiffness matrix of the lower medium
      :type C2: 2D array
      :param theta: incident angle, in degree unit
      :type theta: array-like
      :param azimuth: azimuth angle, in degree unit
      :type azimuth: array-like

      :returns: *float or array-like* -- PP reflectivities















      ..
          !! processed by numpydoc !!

   .. py:method:: Aki_Richards(theta, vp1, vp2, vs1, vs2, den1, den2)
      :staticmethod:

      
      Aki-Richard approximation to PP reflectivity.

      :param theta: incident angle, degree
      :type theta: float or array-like
      :param vp1: P wave velocity of layer 1, m/s
      :type vp1: float
      :param vp2: P wave velocity of layer 2, m/s
      :type vp2: float
      :param vs1: S wave velocity of layer 1, m/s
      :type vs1: float
      :param vs2: S wave velocity of layer 2, m/s
      :type vs2: float
      :param den1: density of layer 1, kg/m3
      :type den1: float
      :param den2: density of layer 2, kg/m3
      :type den2: float

      :returns: *float or array-like* -- R_pp: P wave reflectivity
                R_ps: PS reflectivity
                Rpp0: intercept
                gradient















      ..
          !! processed by numpydoc !!

   .. py:method:: zoeppritz(vp1, vs1, rho1, vp2, vs2, rho2, theta)
      :staticmethod:

      
      Reflection & Transmission coefficients calculated using full Zoeppritz equations.

      :param vp1: P wave velocity of layer 1, m/s
      :type vp1: float
      :param vs1: S wave velocity of layer 1, m/s
      :type vs1: float
      :param rho1: density of layer 1, kg/m3
      :type rho1: float
      :param vp2: P wave velocity of layer 2, m/s
      :type vp2: float
      :param vs2: S wave velocity of layer 2, m/s
      :type vs2: float
      :param rho2: density of layer 2, kg/m3
      :type rho2: float
      :param theta: incident angle, degree
      :type theta: float or array-like

      :returns: *float or array-like* -- Rpp,Rps: PP and PS reflectivity















      ..
          !! processed by numpydoc !!

   .. py:method:: AVO_abe(vp1, vs1, d1, vp2, vs2, d2)
      :staticmethod:

      
      Different approximations AVO terms

      :param vp1: P wave velocity of layer 1, m/s
      :type vp1: float or array-like
      :param vs1: S wave velocity of layer 1, m/s
      :type vs1: float or array-like
      :param d1: density of layer 1, kg/m3
      :type d1: float or array-like
      :param vp2: P wave velocity of layer 2, m/s
      :type vp2: float or array-like
      :param vs2: S wave velocity of layer 2, m/s
      :type vs2: float or array-like
      :param d2: density of layer 2, kg/m3
      :type d2: float or array-like

      :returns: *float or array-like* -- different linear AVO approximations















      ..
          !! processed by numpydoc !!

   .. py:method:: EI_ref(Vp, Vs, rho, theta, SP, norm=True)
      :staticmethod:

      
      Compute elastic impedance of an isotropic, flat-layered Earth

      :param vp1: P wave velocity of layer 1, m/s
      :type vp1: float or array-like
      :param vs1: S wave velocity of layer 1, m/s
      :type vs1: float or array-like
      :param d1: density of layer 1, kg/m3
      :type d1: float or array-like
      :param Vp: P wave velocity
      :type Vp: float or array-like
      :param Vs: S wave velocity
      :type Vs: float or array-like
      :param rho: density
      :type rho: float or array-like
      :param theta: incident angles
      :type theta: array-like
      :param SP: constant ratio of Vs to Vp, can be taken as the average of input Vs/Vp, i.e. SP= VS.mean()/VP.mean()
      :type SP: float
      :param norm: If True: normalized input velocities and density such that the units and dimension match with acoustic impedance. Defaults to True.
      :type norm: bool, optional

      :returns: *float or array-like* -- EI_pp: elastic impedance for PP reflection
                EI_svp: elastic impedance for P-SV reflection
                EI_psv: elastic impedance for SV-P reflection
                EI_svsv: elastic impedance for SV-SV reflection
                EI_shsh: elastic impedance for SH-SH reflection















      ..
          !! processed by numpydoc !!

   .. py:method:: AVO_ortho(a1, b1, e11, d11, e12, d12, g1, rho1, a2, b2, e21, d21, e22, d22, g2, rho2, the)
      :staticmethod:

      
      calculates the reflectivity in the symmetry plane for interfaces between 2 orthorhombic media, refactered from srb toolbox written by Diana Sava.
      :param a1: P-wave vertical velocities of upper medium (1)
      :type a1: float or array-like
      :param b1: S-wave vertical velocities of upper medium (1)
      :type b1: float or array-like
      :param e11: epsilon in the two symmetry planes of the orthorhombic medium for the upper medium (first index indicates the upper medium (1), second index indicates the plane of symmetry (1 - plane perpendicular to x, 2 - plane perpendicular to y);
      :type e11: float or array-like
      :param d11: delta in the two symmetry planes of the orthorhombic medium for the upper medium
      :type d11: float or array-like
      :param e12: epsilon in the two symmetry planes of the orthorhombic medium for the upper medium
      :type e12: float or array-like
      :param d12: delta in the two symmetry planes of the orthorhombic medium for the upper medium
      :type d12: float or array-like
      :param g1: vertical shear wave splitting parameter for the upper medium (1)
      :type g1: float or array-like
      :param rho1: density of the upper medium
      :type rho1: float or array-like
      :param a2: P-wave vertical velocities of lower medium (2)
      :type a2: float or array-like
      :param b2: S-wave vertical velocities of lower medium (2)
      :type b2: float or array-like
      :param e21: epsilon in the two symmetry planes of the orthorhombic medium for the lower medium
      :type e21: float or array-like
      :param d21: delta in the two symmetry planes of the orthorhombic medium for the lower medium
      :type d21: float or array-like
      :param e22: epsilon in the two symmetry planes of the orthorhombic medium for the lower medium
      :type e22: float or array-like
      :param d22: delta in the two symmetry planes of the orthorhombic medium for the lower medium
      :type d22: float or array-like
      :param g2: vertical shear wave splitting parameter for the upper medium (2)
      :type g2: float or array-like
      :param rho2: density of the lower medium
      :type rho2: float or array-like
      :param the: incident angle
      :type the: float or array-like

      :returns: *array-like* -- Rxy: PP reflectivity as a function of angle of incidence in xz plane (13).
                Ryz: PP reflectivity as a function of angle of incidence in yz plane (23)















      ..
          !! processed by numpydoc !!


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


.. py:class:: GM

   
   Contact based granular medium models and extensions.
















   ..
       !! processed by numpydoc !!
   .. py:method:: ThomasStieber(phi_sand, phi_sh, vsh)
      :staticmethod:

      
      Thomas-Stieber porosity model for sand-shale system.

      :param phi_sand: clean sand porosity
      :type phi_sand: float
      :param phi_sh: shale porosity
      :type phi_sh: float
      :param vsh: volume faction of shale in the mixture
      :type vsh: float or array-like

      :returns: *float or array-like* -- phi_ABC,phi_AC (frac): porosity line as shown in Fig 5.3.2 in (Mavko,2020)















      ..
          !! processed by numpydoc !!

   .. py:method:: silty_shale(C, Kq, Gq, Ksh, Gsh)
      :staticmethod:

      
      Dvorkin–Gutierrez silty shale model: model the elastic moduli of decreasing clay content for shale.

      :param C: volume fraction of clay
      :type C: float or array-like
      :param Kq: bulk modulus of silt grains
      :type Kq: float
      :param Gq: shear modulus of silt grains
      :type Gq: float
      :param Ksh: saturated bulk modulus of pure shale
      :type Ksh: float
      :param Gsh: saturated shear modulus of pure shale, * Ksh and Gsh could be derived from well-log measurements of VP, VS and density in a pure shale zone.
      :type Gsh: float

      :returns: *float or array-like* -- K_sat, G_sat: elastic moduli of the saturated silty shale.















      ..
          !! processed by numpydoc !!

   .. py:method:: shaly_sand(phis, C, Kss, Gss, Kcc, Gcc)
      :staticmethod:

      
      Modeling elastic moduli for sand with increasing clay content using LHS bound rather than using Gassmann relation.

      :param phis: critical porosity of sand composite
      :type phis: float
      :param C: clay content
      :type C: float or array-like
      :param Kss: saturated bulk moduli for clean sandstone using e.g. HM
      :type Kss: float
      :param Gss: saturated shear moduli for clean sandstone using e.g. HM
      :type Gss: float
      :param Kcc: saturated bulk moduli calculated from the sandy shale model at critical clay content using silty shale model
      :type Kcc: float
      :param Gcc: saturated shear moduli calculated from the sandy shale model at critical clay content using silty shale model
      :type Gcc: float

      :returns: *float or array-like* -- K_sat,G_sat: saturated rock moduli of the shaly sand















      ..
          !! processed by numpydoc !!

   .. py:method:: contactcement(K0, G0, Kc, Gc, phi, phic, Cn, scheme)
      :staticmethod:

      
      Compute dry elastic moduli of cemented sandstone via Contact cement model by Dvorkin &Nur (1996).

      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param Kc: Bulk modulus of cement
      :type Kc: float
      :param Gc: Shear modulus of cement
      :type Gc: float
      :param phi: Porosity
      :type phi: float or array-like
      :param phic: Critical Porosity
      :type phic: float
      :param Cn: coordination number
      :type Cn: float
      :param scheme:
                     Scheme of cement deposition
                             1=cement deposited at grain contacts
                             2=cement deposited at grain surfaces
      :type scheme: int

      :returns: *_type_* -- K_dry, G_dry (GPa): Effective elastic moduli of dry rock

      .. rubric:: References

      - Dvorkin & Nur, 1996, Geophysics, 61, 1363-1370















      ..
          !! processed by numpydoc !!

   .. py:method:: hertzmindlin(K0, G0, phic, Cn, sigma, f)
      :staticmethod:

      
      Compute effective dry elastic moduli of granular packing under hydrostatic pressure condition via Hertz-Mindlin approach. Reduced shear factor that honours the non-uniform contacts in the granular media is implemented.

      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
                 phic : float
                 Critical Porosity
      :type G0: float
      :param Cn: coordination number
      :type Cn: float
      :param sigma: effective stress
      :type sigma: float or array-like
      :param f: reduced shear factor between 0 and 1
                0=dry pack with inifinitely rough spheres;
                1=dry pack with infinitely smooth spheres
      :type f: float

      :returns: * **K_dry, G_dry** (*float or array-like*) -- effective elastic moduli of dry pack
                * *References* -- ----------
                * *- Rock physics handbook section 5.5.*
                * *- Bachrach, R. and Avseth, P. (2008) Geophysics, 73(6), E197–E209.*















      ..
          !! processed by numpydoc !!

   .. py:method:: softsand(K0, G0, phi, phic, Cn, sigma, f)
      :staticmethod:

      
      Soft-sand (unconsolidated sand) model: model the porosity-sorting effects using the lower Hashin-Shtrikman-Walpole bound. (Also referred to as the 'friable-sand model' in Avseth et al. (2010).

      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param phi: Porosity
                  phic : float
                  Critical Porosity
      :type phi: float or array like
      :param Cn: coordination number
      :type Cn: float
      :param sigma: effective stress
      :type sigma: float or array-like
      :param f: reduced shear factor between 0 and 1
                0=dry pack with inifinitely rough spheres;
                1=dry pack with infinitely smooth spheres
      :type f: float

      :returns: * *float or array-like* -- K_dry, G_dry (GPa): Effective elastic moduli of dry pack
                * *References* -- ----------
                * *- The Uncemented (Soft) Sand Model in Rock physics handbook section 5.5*
                * *- Avseth, P.; Mukerji, T. & Mavko, G. Cambridge university press, 2010*















      ..
          !! processed by numpydoc !!

   .. py:method:: Walton(K0, G0, phic, Cn, sigma, f)
      :staticmethod:

      
      Compute dry rock elastic moduli of sphere packs based on the Walton (1987)' thoery. Reduced shear factor that honours the non-uniform contacts in the granular media is implemented.

      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
                 phic : float
                 Critical Porosity
      :type G0: float
      :param Cn: coordination number
      :type Cn: float
      :param sigma: effective stress
      :type sigma: float or array-like
      :param f: reduced shear factor between 0 and 1
                0=dry pack with inifinitely rough spheres;
                1=dry pack with infinitely smooth spheres
      :type f: float

      :returns: * *float or array-like* -- K_w, G_w: Effective elastic moduli of dry pack
                * *References* -- ----------
                * *- Walton model in Rock physics handbook section 5.5*
                * *- Walton, K., 1987, J. Mech. Phys. Solids, vol.35, p213-226.*
                * *- Bachrach, R. and Avseth, P. (2008) Geophysics, 73(6), E197–E209*















      ..
          !! processed by numpydoc !!

   .. py:method:: johnson(K0, G0, n, phi, epsilon, epsilon_axial, path='together')
      :staticmethod:

      
      effective theory for stress-induced anisotropy in sphere packs. The transversely isotropic strain is considered as a combination of hydrostatic strain and uniaxial strain.

      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param n: coordination number
      :type n: float
      :param phi: porosity
      :type phi: float or array like
      :param epsilon: hydrostatic strain (negative in compression)
      :type epsilon: float or array like
      :param epsilon_axial: uniaxial strain (along 3-axis)
      :type epsilon_axial: float or array like
      :param path: 'together': the hydrostatic and uniaxial strains are applied simultaneously
                   'uni_iso': the uniaxial strain is applied first followed by a hydrostatic strain
                   'iso_uni': the hydrostatic strain is applied first followed by a uniaxial strain by default 'together'
      :type path: str, optional

      :returns: * *array and float* -- C: (matrix): VTI stiffness matrix
                  sigma33: non zero stress tensor component
                  sigma11: non zero stress tensor component, sigma11=sigma22
                * *References* -- ----------
                * *- Norris, A. N., and Johnson, D. L., 1997, ASME Journal of Applied Mechanics, 64, 39-49.*
                * *- Johnson, D.L., Schwartz, L.M., Elata, D., et al., 1998. Transactions ASME, 65, 380–388.*















      ..
          !! processed by numpydoc !!

   .. py:method:: stiffsand(K0, G0, phi, phic, Cn, sigma, f)
      :staticmethod:

      
      Stiff-sand model:  Modified Hashin-Shtrikman upper bound with Hertz-Mindlin end point, counterpart to soft sand model.
      model the porosity-sorting effects using the lower Hashin–Shtrikman–Walpole bound.

      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param phi: Porosity
                  phic : float
                  Critical Porosity
      :type phi: float or array like
      :param Cn: coordination number
      :type Cn: float
      :param sigma: effective stress
      :type sigma: float or array-like
      :param f: reduced shear factor between 0 and 1
                0=dry pack with inifinitely rough spheres;
                1=dry pack with infinitely smooth spheres
      :type f: float

      :returns: *float or array-like* -- K_dry, G_dry (GPa): Effective elastic moduli of dry pack















      ..
          !! processed by numpydoc !!

   .. py:method:: constantcement(phi_b, K0, G0, Kc, Gc, phi, phic, Cn, scheme)
      :staticmethod:

      
      Constant cement (constant depth) model according to Avseth (2000)

      :param phi_b: adjusted high porosity end memeber
      :type phi_b: _type_
      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param Kc: Bulk modulus of cement
      :type Kc: float
      :param Gc: Shear modulus of cement
      :type Gc: float
      :param phi: Porosity
      :type phi: float or array-like
      :param phic: Critical Porosity
      :type phic: float
      :param Cn: coordination number
      :type Cn: float
      :param scheme:
                     Scheme of cement deposition
                             1=cement deposited at grain contacts
                             2=cement deposited at grain surfaces
      :type scheme: int

      :returns: * *float or array-like* -- K_dry, G_dry (GPa): Effective elastic moduli of dry rock
                * *References* -- ----------
                * *- Avseth, P.; Dvorkin, J.; Mavko, G. & Rykkje, J. Geophysical Research Letters, Wiley Online Library, 2000, 27, 2761-2764*















      ..
          !! processed by numpydoc !!

   .. py:method:: MUHS(K0, G0, Kc, Gc, phi, phi_b, phic, Cn, scheme)
      :staticmethod:

      
      Increasing cement model: Modified Hashin-Strikmann upper bound blend with contact cement model. For elastically stiff sandstone modelling.

      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param Kc: Bulk modulus of cement
      :type Kc: float
      :param Gc: Shear modulus of cement
      :type Gc: float
      :param phi: Porosity
      :type phi: float or array-like
      :param phi_b: adjusted high porosity end memeber
      :type phi_b: _type_
      :param phic: Critical Porosity
      :type phic: float
      :param Cn: coordination number
      :type Cn: float
      :param scheme:
                     Scheme of cement deposition
                             1=cement deposited at grain contacts
                             2=cement deposited at grain surfaces
      :type scheme: int

      :returns: * *float or array-like* -- K_dry, G_dry (GPa): Effective elastic moduli of dry rock
                * *References* -- ----------
                * *- Avseth, P.; Mukerji, T. & Mavko, G. Cambridge university press, 2010*















      ..
          !! processed by numpydoc !!

   .. py:method:: Digby(K0, G0, phi, Cn, sigma, a_R)
      :staticmethod:

      
      Compute Keff and Geff using Digby's model

      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param phi: Porosity
      :type phi: float
      :param Cn: coordination number
      :type Cn: float
      :param sigma: stress
      :type sigma: float or array-like
      :param a_R: a_R (unitless): ratio of the radius of the initially bonded area to the grain radius
      :type a_R: float

      :returns: * *float or array-like* -- Keff, Geff (Gpa): effective medium stiffness
                * *References* -- ----------
                * *- Digby, P.J., 1981. Journal of Applied Mechanics, 48, 803–808.*















      ..
          !! processed by numpydoc !!

   .. py:method:: pcm(f, sigma, K0, G0, phi, phic, v_cem, v_ci, Kc, Gc, Cn, mode, scheme, f_)
      :staticmethod:

      
      Computes effective elastic moduli of patchy cemented sandstone according to Avseth (2016).

      :param f: volume fraction of cemented rock in the binary mixture
      :type f: float
      :param sigma: effective stress
      :type sigma: float or array-like
      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param phi: Porosity
      :type phi: float
      :param phic: Critical Porosity
      :type phic: float
      :param v_cem: cement fraction in contact cement model. phi_cem= phic-vcem
      :type v_cem: float
      :param v_ci: cement threshold above which increasing cement model is applied
      :type v_ci: float
      :param Kc: bulk modulus of cement
      :type Kc: float
      :param Gc: shear modulus of cement
      :type Gc: float
      :param Cn: coordination number
      :type Cn: float
      :param mode: 'stiff' or 'soft'. stiffest mixing or softest mixing. Defaults to 'stiff'.
      :type mode: str
      :param scheme: contact cement scheme.
                     1=cement deposited at grain contacts
                     2=cement deposited at grain surfaces
      :type scheme: int
      :param f_: slip factor in HM modelling. Defaults to 0.5.
      :type f_: float

      .. note:: (Avseth,2016): If 10% is chosen as the “critical” cement limit, the  increasing cement model can be used in addition to the contact cement model. (Torset, 2020): with the increasing cement model appended at 4% cement"

      :returns: * *float or array-like* -- K_DRY, G_DRY (GPa): effective elastic moduli of the dry rock
                * *References* -- ----------
                  - Avseth, P.; Skjei, N. & Mavko, G. The Leading Edge, GeoScienceWorld, 2016, 35, 868-87.















      ..
          !! processed by numpydoc !!

   .. py:method:: diluting(k, sigma0, sigma, m)
      :staticmethod:

      
      stress dependent diluting parameter used in varying patchiness cement model.

      :param k: cement crushing factor. k<=1: no cement crumbling; k>1: cement crumbling.
      :type k: float
      :param sigma0: reference stress, e.g. maximum effective stress, stress at which unloading begins.
      :type sigma0: float
      :param sigma: effective stress
      :type sigma: array-like
      :param m: curvature parameter that defines diluting rate.
      :type m: float

      :returns: *array-like* -- stress dependent diluting parameter















      ..
          !! processed by numpydoc !!

   .. py:method:: vpcm(alpha, f, sigma, K0, G0, phi, phic, v_cem, v_ci, Kc, Gc, Cn, scheme, f_)
      :staticmethod:

      
      Compute effective elastic moduli using varying patchiness cement model (VPCM) as proposed by Yu et al. (2023).

      :param alpha: diluting parameters
      :type alpha: float or array-like
      :param f: volume fraction of cemented rock in the binary mixture
      :type f: float
      :param sigma: effective stress
      :type sigma: float or array-like
      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param phi: Porosity
      :type phi: float
      :param phic: Critical Porosity
      :type phic: float
      :param v_cem: cement fraction in contact cement model. phi_cem= phic-vcem
      :type v_cem: float
      :param v_ci: cement threshold above which increasing cement model is applied
      :type v_ci: float
      :param Kc: bulk modulus of cement
      :type Kc: float
      :param Gc: shear modulus of cement
      :type Gc: float
      :param Cn: coordination number
      :type Cn: float
      :param scheme: contact cement scheme.
                     1=cement deposited at grain contacts
                     2=cement deposited at grain surfaces
      :type scheme: int
      :param f_: slip factor in HM modelling.
      :type f_: float
      :param Note: (Avseth,2016): If 10% is chosen as the “critical” cement limit, the increasing cement model can be used in addition to the contact cement model. (Torset, 2020): with the increasing cement model appended at 4% cement"

      :returns: *array-like* -- K_DRY, G_DRY (GPa): effective elastic moduli of the dry rock















      ..
          !! processed by numpydoc !!


.. py:class:: Fluid

   
   Fluid subsitution approaches and models describing velocity dispersion and attenuation due to the fluid effect.
















   ..
       !! processed by numpydoc !!
   .. py:method:: Brie(Kw, Kgas, Sw, e)
      :staticmethod:

      
      Brie empirical fluid mixing law

      :param Kw: bulk modulus of fluid phase
      :type Kw: float
      :param Kgas: bulk modulus of gas phase
      :type Kgas: float
      :param Sw: water saturation
      :type Sw: float or array
      :param e: Brie component
      :type e: int

      :returns: *float or array* -- Kf: effective fluid propertie















      ..
          !! processed by numpydoc !!

   .. py:method:: Biot(Kdry, Gdry, K0, Kfl, rho0, rhofl, eta, phi, kapa, a, alpha, freq)
      :staticmethod:

      
      Compute Biot dispersion and velocity attenuation

      :param Kdry: dry frame bulk modulus
      :type Kdry: float or array-like
      :param Gdry: dry frame shear modulus
      :type Gdry: float or array-like
      :param K0: bulk modulus of mineral material making up rock
      :type K0: float
      :param Kfl: effective bulk modulus of pore fluid
      :type Kfl: float
      :param rho0: grain density
      :type rho0: float
      :param rhofl: pore fluid density
      :type rhofl: float
      :param eta: η is the viscosity of the pore fluid
      :type eta: float
      :param phi: porosity
      :type phi: float
      :param kapa: absolute permeability of the rock
      :type kapa: float
      :param a: pore-size parameter. Stoll (1974) found that values between 1/6 and 1/7 of the mean grain diameter
      :type a: float
      :param alpha: tortuosity parameter, always greater than or equal to 1.
      :type alpha: float
      :param freq: frequency range, e.g 10^-3 to 10^3 Hz
      :type freq: float or array-like

      :returns: *float or array-like* -- Vp_fast, fast P-wave velocities at all frequencies,km/s
                Vp_slow, slow P-wave velocities at all frequencies,km/s
                Vs, S-wave velocities,km/s
                QP1_inv, fast P-wave attenuation
                QP2_inv, slow P-wave attenuation
                Qs_inv, S-wave attenuation















      ..
          !! processed by numpydoc !!

   .. py:method:: Biot_HF(Kdry, Gdry, K0, Kfl, rho0, rhofl, phi, alpha)
      :staticmethod:

      
      Biot high-frequency limiting velocities in the notation of Johnson and Plona (1982)

      :param Kdry: dry frame bulk modulus
      :type Kdry: float or array-like
      :param Gdry: dry frame shear modulus
      :type Gdry: float or array-like
      :param K0: bulk modulus of mineral material making up rock
      :type K0: float
      :param Kfl: effective bulk modulus of pore fluid
      :type Kfl: float
      :param rho0: grain density
      :type rho0: float
      :param rhofl: pore fluid density
      :type rhofl: float
      :param phi: porosity
      :type phi: float
      :param alpha: tortuosity parameter, always greater than or equal to 1.
      :type alpha: float

      :returns: *float or array-like* -- Vp_fast,Vp_slow,Vs:  high-frequency limiting velocities,km/s















      ..
          !! processed by numpydoc !!

   .. py:method:: Geertsma_Smit_HF(Kdry, Gdry, K0, Kfl, rho0, rhofl, phi, alpha)
      :staticmethod:

      
      Approximation of Biot high-frequency limit of the fast P-wave velocity given by Geertsma and Smit (1961), This form predicts velocities that are too high (by about 3%–6%) compared with the actual high-frequency limit.

      :param Kdry: dry frame bulk modulus
      :type Kdry: float or array-like
      :param Gdry: dry frame shear modulus
      :type Gdry: float or array-like
      :param K0: bulk modulus of mineral material making up rock
      :type K0: float
      :param Kfl: effective bulk modulus of pore fluid
      :type Kfl: float
      :param rho0: grain density
      :type rho0: float
      :param rhofl: pore fluid density
      :type rhofl: float
      :param phi: porosity
      :type phi: float
      :param alpha: tortuosity parameter, always greater than or equal to 1.
      :type alpha: float

      :returns: *float or array-like* -- Vp_fast,Vs: high-frequency limiting velocities, km/s















      ..
          !! processed by numpydoc !!

   .. py:method:: Geertsma_Smit_LF(Vp0, Vpinf, freq, phi, rhofl, kapa, eta)
      :staticmethod:

      
      Low and middle-frequency approximations of Biot wave given by Geertsma and Smit (1961). Noticed that mathematically this approximation is valid at moderate-to-low seismic frequencies, i.e. f<fc

      :param Vp0: Biot−Gassmann low-frequency limiting P-wave velocity, km/s or m/s
      :type Vp0: float
      :param Vpinf: Biot highfrequency limiting P-wave velocity, km/s or m/s
      :type Vpinf: float
      :param freq: frequency
      :type freq: float or array-like
      :param phi: porosity
      :type phi: float
      :param rhofl: fluid density
      :type rhofl: float
      :param kapa: absolute permeability of the rock.
      :type kapa: float
      :param eta: viscosity of the pore fluid
      :type eta: float

      :returns: *float or array-like* -- Vp: frequency-dependent P-wave velocity of saturated rock















      ..
          !! processed by numpydoc !!

   .. py:method:: Gassmann(K_dry, G_dry, K_mat, Kf, phi)
      :staticmethod:

      
      Computes saturated elastic moduli of rock via Gassmann equation given dry-rock moduli.

      :param K_dry: dry frame bulk modulus
      :type K_dry: float or array-like
      :param G_dry: dry frame shear modulus
      :type G_dry: float or array-like
      :param K_mat: matrix bulk modulus
      :type K_mat: float
      :param Kf: fluid bulk modulus
      :type Kf: float
      :param phi: porosity
      :type phi: float or array-like

      :returns: *float or array-like* -- K_sat, G_sat: fluid saturated elastic moduli















      ..
          !! processed by numpydoc !!

   .. py:method:: Gassmann_sub(phi, K0, Ksat1, Kfl1, Kfl2)
      :staticmethod:

      
      Fluid subsititution using Gassmann equation, thr rock is initially saturated with a fluid, compute the saturated moduli for tge rock saturated with a different fluid

      :param phi: porosity
      :type phi: float or array-like
      :param K0: mineral modulus
      :type K0: float
      :param Ksat1: original bulk modulus of rock saturated with fluid of bulk modulus Kfl1
      :type Ksat1: float or array-like
      :param Kfl1: original saturant
      :type Kfl1: float
      :param Kfl2: new saturant
      :type Kfl2: float

      :returns: *float or array-like* -- Ksat2: new satuarted bulk modulus of the rock















      ..
          !! processed by numpydoc !!

   .. py:method:: vels(K_dry, G_dry, K0, D0, Kf, Df, phi)
      :staticmethod:

      
      Computes Vp,Vs and densities of saturated rock using Gassmann relations from elastic moduli of rock. See also `Gassmann_vels`.

      :param K_dry: dry frame bulk modulus
      :type K_dry: float
      :param G_dry: dry frame shear modulus
      :type G_dry: float
      :param K0: mineral matrix bulk modulus
      :type K0: float
      :param D0: mineral matrix density
      :type D0: float
      :param Kf: fluid bulk modulus
      :type Kf: float
      :param Df: fluid density in g/cm3
      :type Df: float
      :param phi: porosity
      :type phi: float or array

      :returns: *float or array* -- Vp, Vs, rho















      ..
          !! processed by numpydoc !!

   .. py:method:: Gassmann_vels(Vp1, Vs1, rho1, rhofl1, Kfl1, rhofl2, Kfl2, K0, phi)
      :staticmethod:

      
      Gassmann fluid substituion with velocities

      :param Vp1: saturated P velocity of rock with fluid 1
      :type Vp1: float or array-like
      :param Vs1: saturated S velocity of rock with fluid 1
      :type Vs1: float or array-like
      :param rho1: bulk density of saturated rock with fluid 1
      :type rho1: float
      :param rhofl1: density of fluid 1
      :type rhofl1: float
      :param Kfl1: bulk modulus of fluid 1
      :type Kfl1: float
      :param rhofl2: density of fluid 2
      :type rhofl2: float
      :param Kfl2: bulk modulus of fluid 2
      :type Kfl2: float
      :param K0: mineral bulk modulus
      :type K0: float
      :param phi: porosity
      :type phi: float or array-like

      :returns: *float or array-like* -- Vp2, Vs2: velocities of rock saturated with fluid 2















      ..
          !! processed by numpydoc !!

   .. py:method:: Gassmann_approx(Msat1, M0, Mfl1, phi, Mfl2)
      :staticmethod:

      
      Perform gassmann fluid subsititution using p wave modulus only

      :param Msat1: in situ saturated p wave modulus from well log data
      :type Msat1: float or array-like
      :param M0: p wave modulus of mineral
      :type M0: float
      :param Mfl1: p wave modulus of in situ fluid
      :type Mfl1: float
      :param phi: porosity
      :type phi: float
      :param Mfl2: p wave modulus of new fluid for susbtitution
      :type Mfl2: float

      :returns: *float or array-like* -- Msat2:  p wave modulus of rock fully saturated with new fluid















      ..
          !! processed by numpydoc !!

   .. py:method:: Brown_Korringa_dry2sat(Sdry, K0, G0, Kfl, phi)
      :staticmethod:

      
      Compute fluid saturated compliances from dry compliance for anisotropic rock using Brown and Korringa (1975). See eq. 32 in the paper.

      :param Sdry: comliance matrix of the dry rock
      :type Sdry: 2d array
      :param K0: Isotropic mineral bulk modulus
      :type K0: float
      :param G0: Isotropic mineral shear modulus
      :type G0: float
      :param Kfl: Isotropic fluid bulk modulus
      :type Kfl: float
      :param phi: porosity
      :type phi: float

      :returns: *2d array* -- Ssat (6x6 matrix): Saturated compliance of anisotropic rock















      ..
          !! processed by numpydoc !!

   .. py:method:: Brown_Korringa_sat2dry(Ssat, K0, G0, Kfl, phi)
      :staticmethod:

      
      Compute dry compliance from fluid saturated compliances for arbitrarily anisotropic rock using Brown and Korringa (1975). See eq. 32 in the paper.

      :param Ssat: comliance matrix (6x6) of the saturated rock
      :type Ssat: 2d array
      :param K0: Isotropic mineral bulk modulus
      :type K0: float
      :param G0: Isotropic mineral shear modulus
      :type G0: float
      :param Kfl: Isotropic fluid bulk modulus
      :type Kfl: float
      :param phi: porosity
      :type phi: float

      :returns: *2d array* -- Sdry (6x6 matrix): Dry compliance of anisotropic rock















      ..
          !! processed by numpydoc !!

   .. py:method:: Brown_Korringa_sub(Csat, K0, G0, Kfl1, Kfl2, phi)
      :staticmethod:

      
      Fluid substitution in arbitrarily anisotropic rock using Brown and Korringa (1975). the rock is originally saturated by fluid 1. After fluid subsititution, the rock is finally saturated by fluid 2.

      :param Csat: comliance matrix of the saturated rock
      :type Csat: 6x6 matrix
      :param K0: Isotropic mineral bulk modulus
      :type K0: float
      :param G0: Isotropic mineral shear modulus
      :type G0: float
      :param Kfl1: bulk modulus of the original fluid
      :type Kfl1: float
      :param Kfl2: bulk modulus of the final fluid
      :type Kfl2: float
      :param phi: porosity
      :type phi: float

      :returns: *2d array* -- Csat2, Ssat2 (6x6 matrix): Dry stiffness and compliance matrix of anisotropic rock saturated with new fluid















      ..
          !! processed by numpydoc !!

   .. py:method:: Mavko_Jizba(Vp_hs, Vs_hs, Vpdry, Vsdry, K0, rhodry, rhofl, Kfl, phi)
      :staticmethod:

      
      Predicting the very high-frequency moduli and velocities of saturated rocks from dry rock properties using the Squirt flow model derived by Mavko and Jizba (1991).

      :param Vp_hs: P wave velocity of the dry rock measured at very high effective pressure in the unit of m/s
      :type Vp_hs: float
      :param Vs_hs: S wave velocity of the dry rock  measured at very high effective pressure in the unit of m/s
      :type Vs_hs: float
      :param Vpdry: P wave velocity of the dry rock measured at different effective pressure in the unit of m/s
      :type Vpdry: array
      :param Vsdry: S wave velocity of the dry rock measured at different effective pressure in the unit of m/s
      :type Vsdry: array
      :param K0: mineral bulk moduli
      :type K0: float
      :param rhodry: bulk density of the dry rock
      :type rhodry: float
      :param rhofl: bulk density of the pore fluid
      :type rhofl: float
      :param Kfl: bulk moduli of the pore fluid
      :type Kfl: float
      :param phi: porosity
      :type phi: float

      :returns: *float or array-like* -- Kuf_sat (float):GPa, predicted high frequency bulk moduli of saturated rock
                Guf_sat (array): GPa, predicted high frequency shear moduli of saturated rock at different pressure
                Vp_hf (array): m/s, predicted high frequency P wave velocities of saturated rock
                Vs_hf (array): m/s, predicted high frequency S wave velocities of saturated rock















      ..
          !! processed by numpydoc !!

   .. py:method:: Squirt_anisotropic(Sdry, Sdry_hp)
      :staticmethod:

      
      Predict wet unrelaxed frame compliances at very high frequency from dry frame compliances for transversely isotropic rocks using theoretical formula derived by Mukerji and Mavko, (1994)

      :param Sdry: dry rock compliances [S11 S12 S13 S33 S44]
      :type Sdry: list or array
      :param Sdry_hp: dry rock compliances at very high effective stress [S11 S12 S13 S33 S44]
      :type Sdry_hp: array

      :returns: *array* -- The wet-frame compliances [S11 S12 S13 S33 S44]















      ..
          !! processed by numpydoc !!

   .. py:method:: White_Dutta_Ode(Kdry, Gdry, K0, phi, rho0, rhofl1, rhofl2, Kfl1, Kfl2, eta1, eta2, kapa, a, sg, freq)
      :staticmethod:

      
      Dispersion and Attenuation of partial saturation using White and Dutta–Odé Model.

      :param Kdry: bulk modulus of the dry rock
      :type Kdry: float
      :param Gdry: shear modulus of the dry rock
      :type Gdry: float
      :param K0: Isotropic mineral bulk modulus
      :type K0: float
      :param phi: porosity
      :type phi: float
      :param rho0: mineral density
      :type rho0: float
      :param rhofl1: density of the fluid opcupying the central sphere
      :type rhofl1: float
      :param rhofl2: density of the fluid opcupying the outer sphere
      :type rhofl2: float
      :param Kfl1: bulk modulus of the fluid opcupying the central sphere
      :type Kfl1: float
      :param Kfl2: bulk modulus of the fluid opcupying the outer sphere
      :type Kfl2: float
      :param eta1: viscousity of the fluid opcupying the central sphere
      :type eta1: float
      :param eta2: viscousity of the fluid opcupying the outer sphere
      :type eta2: float
      :param kapa: absolute permeability of the rock
      :type kapa: float
      :param a: radius of central sphere , sg=a3/b3
      :type a: float
      :param sg: saturation of fluid opcupying the central sphere
      :type sg: float
      :param freq: frequencies
      :type freq: float or array-like

      :returns: *float, array-like* -- Vp (m/s): P wave velocity km/s
                a_w: attenuation coefficient
                K_star: complex bulk modulus















      ..
          !! processed by numpydoc !!


.. py:class:: EM

   
   classical bounds and inclusion models.
















   ..
       !! processed by numpydoc !!
   .. py:method:: VRH(volumes, M)
      :staticmethod:

      
      Computes Voigt, Reuss, and Hill Average Moduli Estimate.

      :param volumes: volumetric fractions of N phases
      :type volumes: list or array-like
      :param M: elastic modulus of the N phase.
      :type M: list or array-like

      :returns: *float* -- M_v: Voigt average
                M_r: Reuss average
                M_0: Hill average















      ..
          !! processed by numpydoc !!

   .. py:method:: cripor(K0, G0, phi, phic)
      :staticmethod:

      
      Critical porosity model according to Nur’s modified Voigt average.

      :param K0: mineral bulk modulus in GPa
      :type K0: float or array-like
      :param G0: mineral shear modulus in GPa
      :type G0: float or array-like
      :param phi: porosity in frac
      :type phi: float or array-like
      :param phic: critical porosity in frac
      :type phic: float

      :returns: *float or array-like* -- K_dry,G_dry (GPa): dry elastic moduli of the framework















      ..
          !! processed by numpydoc !!

   .. py:method:: cripor_reuss(M0, Mf, phic, den=False)
      :staticmethod:

      
      In the suspension domain, the effective bulk and shear moduli of the rock can be estimated by using the Reuss (isostress) average.

      :param M0: The solid phase modulus or density
      :type M0: float
      :param Mf: The pore filled phase modulus or density
      :type Mf: float
      :param phic: critical porosity
      :type phic: float
      :param den: If False: compute the reuss average for effective modulus of two mixing phases. If true, compute avearge density using mass balance, which corresponds to voigt average. Defaults to False.
      :type den: bool, optional

      :returns: *float or array-like* -- M (GPa/g.cc): average modulus or average density

      .. rubric:: References

      - Section 7.1 Rock physics handbook 2nd edition















      ..
          !! processed by numpydoc !!

   .. py:method:: HS(f, K1, K2, G1, G2, bound='upper')
      :staticmethod:

      
      Compute effective moduli of two-phase composite using hashin-strikmann bounds.

      :param f: 0-1, volume fraction of stiff material
      :type f: float
      :param K1: bulk modulus of stiff phase
      :type K1: float or array-like
      :param K2: bulk modulus of soft phase
      :type K2: float or array-like
      :param G1: shear modulus of stiff phase
      :type G1: float or array-like
      :param G2: shear modulus of soft phase
      :type G2: float or array-like
      :param bound: upper bound or lower bound. Defaults to 'upper'.
      :type bound: str, optional

      :returns: *float or array-like* -- K, G (GPa): effective moduli of two-phase composite















      ..
          !! processed by numpydoc !!

   .. py:method:: Eshelby_Cheng(K, G, phi, alpha, Kf, mat=False)
      :staticmethod:

      
      Compute the effective anisotropic moduli of a cracked isotropic rock with single set fracture using Eshelby–Cheng Model.

      :param K: bulk modulus of the isotropic matrix GPa
      :type K: float
      :param G: shear modulus of the isotropic matrix GPa
      :type G: float
      :param phi: (crack) porosity
      :type phi: float
      :param alpha: aspect ratio of crack
      :type alpha: float
      :param Kf: bulk modulus of the fluid. For dry cracks use fluid bulk modulus 0
      :type Kf: float
      :param mat: If true: the output is in matrix form, otherwise  is numpy array. Defaults to False.
      :type mat: bool, optional

      :returns: *1d or 2d array* -- C_eff: effective moduli of cracked, transversely isotropic rocks

      .. rubric:: References

      - section 4.14 in The Rock Physics Handbook















      ..
          !! processed by numpydoc !!

   .. py:method:: hudson(K, G, Ki, Gi, alpha, crd, order=1, axis=3)
      :staticmethod:

      
      Hudson’s effective crack model assuming weak inclusion for media with single crack set with all normals aligned along 1 or 3-axis. First and Second order corrections are both implemented. Notice that the second order correction has limitation. See Cheng (1993).

      :param K: bulk modulus of isotropic background
      :type K: float
      :param G: shear modulus of isotropic background
      :type G: float
      :param Ki: bulk modulus of the inclusion material. For dry cracks: Ki=0
      :type Ki: float
      :param Gi: shear modulus of the inclusion material
      :type Gi: float
      :param alpha: crack aspect ratio
      :type alpha: float
      :param crd: crack density
      :type crd: float
      :param order:
                    approximation order.
                        1: Hudson's model with first order correction.
                        2: Hudson's model with first order correction.
                        Defaults to 1.
      :type order: int, optional
      :param axis:
                   axis of symmetry.
                       1: crack normals aligned along 1-axis, output HTI
                       3: crack normals aligned along 3-axis, output VTI
                       Defaults to 3
      :type axis: int, optional

      :returns: *2d array* -- C_eff: effective moduli in 6x6 matrix form.















      ..
          !! processed by numpydoc !!

   .. py:method:: hudson_rand(K, G, Ki, Gi, alpha, crd)
      :staticmethod:

      
      Hudson's crack model of a material containing randomly oriented inclusions. The model results agree with the consistent results of Budiansky and O’Connell (1976).

      :param K: bulk modulus of isotropic background
      :type K: float or array-like
      :param G: shear modulus of isotropic background
      :type G: float or array-like
      :param Ki: bulk modulus of the inclusion material. For dry cracks: Ki=0
      :type Ki: float
      :param Gi: shear modulus of the inclusion material, for fluid, Gi=0
      :type Gi: float
      :param alpha: crack aspect ratio
      :type alpha: float
      :param crd: crack density
      :type crd: float

      :returns: *float or array-like* -- K_eff, G_eff (GPa): effective moduli of the medium with randomly oriented inclusions















      ..
          !! processed by numpydoc !!

   .. py:method:: hudson_ortho(K, G, Ki, Gi, alpha, crd)
      :staticmethod:

      
      Hudson’s first order effective crack model assuming weak inclusion for media with three crack sets with normals aligned along 1 2, and 3-axis respectively.  Model is valid for small crack density and aspect ratios.

      :param K: bulk modulus of isotropic background
      :type K: float
      :param G: shear modulus of isotropic background
      :type G: float
      :param Ki: bulk modulus of the inclusion material. For dry cracks: Ki=0
      :type Ki: float
      :param Gi: shear modulus of the inclusion material, for fluid, Gi=0
      :type Gi: float
      :param alpha: [alpha1, alpha2,alpha3] aspect ratios of  three crack sets
      :type alpha: nd array with size 3
      :param crd: [crd1, crd2, crd3] crack densities of three crack sets
      :type crd: nd array with size 3

      :returns: *2d array* -- C_eff: effective moduli in 6x6 matrix form.















      ..
          !! processed by numpydoc !!

   .. py:method:: hudson_cone(K, G, Ki, Gi, alpha, crd, theta)
      :staticmethod:

      
      Hudson’s first order effective crack model assuming weak inclusion for media with crack normals randomly distributed at a fixed angle from the TI symmetry axis 3 forming a cone;

      :param K: bulk modulus of isotropic background
      :type K: float
      :param G: shear modulus of isotropic background
      :type G: float
      :param Ki: bulk modulus of the inclusion material. For dry cracks: Ki=0
      :type Ki: float
      :param Gi: shear modulus of the inclusion material, for fluid, Gi=0
      :type Gi: float
      :param alpha: aspect ratios of crack sets
      :type alpha: float
      :param crd: total crack density
      :type crd: float
      :param theta: the fixed angle between the crack normam and the symmetry axis x3. degree unit.
      :type theta: float

      :returns: *2d array* -- C_eff: effective moduli of TI medium in 6x6 matrix form.















      ..
          !! processed by numpydoc !!

   .. py:method:: Berryman_sc(K, G, X, Alpha)
      :staticmethod:

      
      Effective elastic moduli for multi-component composite using Berryman's Consistent (Coherent Potential Approximation) method.See also: PQ_vectorize, Berryman_func

      :param K: 1d array of bulk moduli of N constituent phases, [K1,K2,...Kn]
      :type K: array-like
      :param G: 1d array of shear moduli of N constituent phases, [G1,G2,...Gn]
      :type G: array-like
      :param X: 1d array of volume fractions of N constituent phases, [x1,...xn], Sum(X) = 1.
      :type X: array-like
      :param Alpha: aspect ratios of N constituent phases. Note that α <1 for oblate spheroids and α > 1 for prolate spheroids, α = 1 for spherical pores,[α1,α2...αn]
      :type Alpha: array-like

      :returns: *array-like* -- K_sc,G_sc: Effective bulk and shear moduli of the composite















      ..
          !! processed by numpydoc !!

   .. py:method:: PQ_vectorize(Km, Gm, Ki, Gi, alpha)
      :staticmethod:

      
      compute geometric strain concentration factors P and Q for prolate and oblate spheroids according to Berymann (1980).See also: Berryman_sc, Berryman_func

      :param Km: bulk modulus of matrix phase. For Berryman SC       approach, this corresponds to the effective moduli of the composite.
      :type Km: float
      :param Gm: shear modulus of matrix phase. For Berryman SC approach, this corresponds to the effective moduli of the composite.
      :type Gm: float
      :param Ki: 1d array of bulk moduli of N constituent phases, [K1,K2,...Kn]
      :type Ki: array-like
      :param Gi: 1d array of shear moduli of N constituent phases, [G1,G2,...Gn]
      :type Gi: array-like
      :param alpha: aspect ratios of N constituent phases. Note that α <1 for oblate spheroids and α > 1 for prolate spheroids, α = 1 for spherical pores,[α1,α2...αn]
      :type alpha: array-like

      :returns: *array-like* -- P,Q (array): geometric strain concentration factors, [P1,,,Pn],[Q1,,,Qn]















      ..
          !! processed by numpydoc !!

   .. py:method:: Berryman_func(params, K, G, X, Alpha)
      :staticmethod:

      
      Form the system of equastions to solve. See 4.11.14 and 4.11.15 in Rock physics handbook 2020. See also: Berryman_sc

      :param params: Parameters to solve, K_sc, G_sc
      :param K: 1d array of bulk moduli of N constituent phases, [K1,K2,...Kn]
      :type K: array
      :param G: 1d array of shear moduli of N constituent phases, [G1,G2,...Gn]
      :type G: array
      :param X: 1d array of volume fractions of N constituent phases, [x1,...xn]
      :type X: array
      :param Alpha: aspect ratios of N constituent phases. Note that α <1 for oblate spheroids and α > 1 for prolate spheroids, α = 1 for spherical pores,[α1,α2...αn]
      :type Alpha: array

      :returns: *equation* -- Eqs to be solved















      ..
          !! processed by numpydoc !!

   .. py:method:: Swiss_cheese(Ks, Gs, phi)
      :staticmethod:

      
      Compute effective elastic moduli via "Swiss cheese" model with spherical pores. "Swiss cheese" model assumes a dilute distribution of spherical inclusions embedded in an * *unbounded* * homogenous solid.  It takes the "noninteracting assumption" in which all cavities (pores) are independent so that their contributions can be added.

      :param Ks: Bulk modulus of matrix in GPa
      :type Ks: float or array-like
      :param Gs: Shear modulus of matrix in GPa
      :type Gs: float or array-like
      :param phi: porosity
      :type phi: float or array-like

      :returns: *float or array-like* -- Kdry,Gdry (GPa): effective elastic moduli















      ..
          !! processed by numpydoc !!

   .. py:method:: SC(phi, Ks, Gs, iter_n)
      :staticmethod:

      
      Self-Consistent(SC) model with spherical pores considering the critical porosity and the interaction effect between inclusions.

      :param phi: porosity in frac, note that phi.shape== Ks.shape
      :type phi: float or array-like
      :param Ks: bulk modulus of matrix phase in GPa
      :type Ks: float
      :param Gs: shear modulus of matrix phase in GPa
      :type Gs: float
      :param iter_n: iterations, necessary iterations increases as f increases.
      :type iter_n: int

      :returns: *float or array-like* -- K_eff,G_eff (GPa): effective elastic moduli















      ..
          !! processed by numpydoc !!

   .. py:method:: Dilute_crack(Ks, Gs, cd)
      :staticmethod:

      
      The non-iteracting randomly oriented crack model.

      :param Ks: bulk modulus of uncracked medium in GPa
      :type Ks: float
      :param Gs: shear modulus of uncracked medium in GPa
      :type Gs: float
      :param cd: crack density
      :type cd: float or array-like

      :returns: *float or array-like* -- K_eff,G_eff (GPa): effective elastic moduli















      ..
          !! processed by numpydoc !!

   .. py:method:: OConnell_Budiansky(K0, G0, crd)
      :staticmethod:

      
      O’Connell and Budiansky (1974) presented equations for effective bulk and shear moduli of a cracked medium with randomly oriented dry penny-shaped cracks (in the limiting case when the aspect ratio α goes to 0)

      :param K0: bulk modulus of background medium
      :type K0: float
      :param G0: shear modulus of background medium
      :type G0: float
      :param crd: crack density
      :type crd: float or array-like

      :returns: *float or array-like* -- K_dry,G_dry: dry elastic moduli of cracked medium















      ..
          !! processed by numpydoc !!

   .. py:method:: OConnell_Budiansky_fl(K0, G0, Kfl, crd, alpha)
      :staticmethod:

      
      Saturated effective elastic moduli using the O’Connell and Budiansky Consistent (SC) formulations under the constraints of small aspect ratio cracks with soft-fluid saturation.

      :param K0: bulk modulus of background medium
      :type K0: float
      :param G0: shear modulus of background medium
      :type G0: float
      :param Kfl: bulk modulus of soft fluid inclusion, e.g gas
      :type Kfl: float
      :param crd: crack density
      :type crd: float
      :param alpha: aspect ratio
      :type alpha: float

      :returns: * *float* -- K_sat,G_sat: elastic moduli of cracked background fully saturated by soft fluid.
                * *References* -- ----------
                * *- O’Connell and Budiansky, (1974)*















      ..
          !! processed by numpydoc !!

   .. py:method:: OC_R_funcs(params, crd, nu_0, w)
      :staticmethod:

      
      Form the system of equastions to solve. Given crack density and w, solve for the D and nu_eff simulaneously using equations 23 and 25 in O’Connell and Budiansky, (1974)

      :param params: Parameters to solve
      :param crd: crack density
      :type crd: float
      :param nu_0: Poisson's ratio of background medium
      :type nu_0: float
      :param w: softness indicator of fluid filled crack, w=Kfl/alpha/K0, soft fluid saturation is w is the order of 1
      :type w: float

      :returns: *equation* -- eqs to be solved















      ..
          !! processed by numpydoc !!

   .. py:method:: PQ(Km, Gm, Ki, Gi, alpha)
      :staticmethod:

      
      compute geometric strain concentration factors P and Q for prolate and oblate spheroids according to Berymann (1980). See also PQ_vectorize

      :param Km: Bulk modulus of matrix phase
      :type Km: float
      :param Gm: Shear modulus of matrix phase
      :type Gm: float
      :param Ki: Bulk modulus of inclusion phase
      :type Ki: float
      :param Gi: Shear modulus of inclusion phase
      :type Gi: float
      :param alpha: aspect ratio of the inclusion. Note that α <1 for oblate spheroids and α > 1 for prolate spheroids
      :type alpha: float

      :returns: *float* -- P,Q (unitless): geometric strain concentration factors















      ..
          !! processed by numpydoc !!

   .. py:method:: DEM(y, t, params)
      :staticmethod:

      
      ODE solver tutorial: https://physics.nyu.edu/pine/pymanual/html/chap9/chap9_scipy.html.
















      ..
          !! processed by numpydoc !!

   .. py:method:: Berryman_DEM(Km, Gm, Ki, Gi, alpha, phi)
      :staticmethod:

      
      Compute elastic moduli of two-phase composites by incrementally adding inclusions of one phase (phase 2) to the matrix phase using Berryman DEM theory

      :param Km: host mineral bulk modulus
      :type Km: float
      :param Gm: host mineral shear modulus
      :type Gm: float
      :param Ki: bulk modulus of inclusion
      :type Ki: float
      :param Gi: shear modulus of inclusion
      :type Gi: float
      :param alpha: aspect ratio of the inclusion phase
      :type alpha: float
      :param phi: desired fraction occupied by the inclusion
      :type phi: float















      ..
          !! processed by numpydoc !!

   .. py:method:: SC_dilute(Km, Gm, Ki, Gi, f, mode)
      :staticmethod:

      
      Elastic solids with elastic micro-inclusions. Random distribution of dilute spherical micro-inclusions in a two phase composite.

      :param Km: bulk modulus of matrix
      :type Km: float
      :param Gm: shear modulus of matrix
      :type Gm: float
      :param Ki: bulk modulus of inclusion
      :type Ki: float
      :param Gi: shear modulus of inclusion
      :type Gi: float
      :param f: volume fraction of inclusion phases
      :type f: float or array
      :param mode: 'stress' if macro stress is prescribed. 'strain' if macro strain is prescribed.
      :type mode: string

      .. rubric:: References

      S. Nemat-Nasser and M. Hori (book) : Micromechanics: Overall Properties of Heterogeneous Materials. Sec 8

      :returns: *float or array* -- K, G: effective moduli of the composite















      ..
          !! processed by numpydoc !!

   .. py:method:: SC_flex(f, iter_n, Km, Ki, Gm, Gi)
      :staticmethod:

      
      iteratively solving self consistent model for a two phase compposite consisting random distribution of spherical inclusion, not limited to pore.

      :param f: volumetric fraction, f.shape== Km.shape
      :type f: float or array
      :param iter_n: iterations, necessary iterations increases as f increases.
      :type iter_n: int
      :param Km: bulk modulus of matrix phase
      :type Km: float
      :param Ki: bulk modulus of inclusion phase
      :type Ki: float
      :param Gm: shear modulus of matrix phase
      :type Gm: float
      :param Gi: shear modulus of inclusion phase
      :type Gi: float
      :param Reference:
      :param ---------:
      :param S. Nemat-Nasser and M. Hori (book):
      :type S. Nemat-Nasser and M. Hori (book): Micromechanics: Overall Properties of Heterogeneous Materials. Sec 8

      :returns: *float or array* -- K_eff, G_eff (GPa): Effective elastic moduli















      ..
          !! processed by numpydoc !!

   .. py:method:: MT_average(f, Kmat, Gmat, K1, G1, K2, G2)
      :staticmethod:

      
      Compute Two-phase composite without matrix using modified Mori-Takana Scheme according to  Iwakuma 2003, one of the inhomogeneities must be considered as a matrix in the limiting model.

      :param f: Volume fraction of matrix/inhomogeneity 1. f1=1-f2, (1-f) can be regarded as pseudo crack density.
      :type f: float or array
      :param Kmat: Bulk modulus of matrix/inhomogeneity 1
      :type Kmat: float
      :param Gmat: shear modulus of  matrix/ inhomogeneity 1
      :type Gmat: float
      :param K1: Bulk modulus of inhomogeneity 1
      :type K1: float
      :param G1: shear modulus of inhomogeneity 1
      :type G1: float
      :param K2: Bulk modulus of inhomogeneity 2
      :type K2: float
      :param G2: shear modulus of inhomogeneity 2
      :type G2: float

      :returns: *float or array* -- K_ave, G_ave [GPa]: MT average bulk and shear modulus















      ..
          !! processed by numpydoc !!


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


.. py:class:: QI(Vp, phi, **kwargs)

   Useful functionalities for quantitative intepretation and relevant tasks.


   initialize the parameters for various QI plots. e.g., phi,den,Vsh,eff_stress, TVD can be given depending on required the input parameters to the plot funtion.

   :param Vp: p wave velocity
   :type Vp: array
   :param phi: porosity
   :type phi: array















   ..
       !! processed by numpydoc !!

   .. py:method:: matrix_modulus(vsh, phi_c, phi_0, M_sh, M_g, M_c)
      :staticmethod:

      
      Calculate the modulus of rock matrix as a function of porosity variation caused by cementation, Note that the rock matrix contains everything excluding pore space.

      :param vsh: bulk volume fraction of shale, can be derived from well log.
      :type vsh: float
      :param phi_c: critical porosity
      :type phi_c: float
      :param phi_0: static porosity during cementation ranging between 0.4 to 0 should be phi when phi is an array of porosity
      :type phi_0: float or array
      :param M_sh: modulus of shale
      :type M_sh: float
      :param M_g: modulus of grain material
      :type M_g: float
      :param M_c: modulus of cement
      :type M_c: float

      :returns: *float or array* -- M_mat: updated matrix modulus















      ..
          !! processed by numpydoc !!

   .. py:method:: den_matrix(vsh, phi_c, phi_0, den_sh, den_g, den_c)
      :staticmethod:

      
      Calculate the matrix density as a function of porosity variation caused by cementation.

      :param vsh: bulk volume fraction of shale, can be derived from well log.
      :type vsh: float
      :param phi_c: critical porosity
      :type phi_c: float
      :param phi_0: static porosity during cementation ranging between 0.4 to 0 should be phi when phi is an array of porosity
      :type phi_0: float or array
      :param den_sh: density of the clay
      :type den_sh: float
      :param den_g: density of the grain
      :type den_g: float
      :param den_c: density of the cement
      :type den_c: float

      :returns: *float or array* -- den_mat: updated matrix density















      ..
          !! processed by numpydoc !!

   .. py:method:: screening(Dqz, Kqz, Gqz, Dsh, Ksh, Gsh, Dc, Kc, Gc, Db, Kb, phib_p, phi_c, sigma, vsh, scheme, f, Cn)
      :staticmethod:

      
      compute elastic bounds used for rock physics screening, the lower bound is computed using friable sand model, and the upper bound is contact cement model blend with increasing cement model.

      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param Dqz: Density of the grain. not limited to quartz grain
      :type Dqz: float
      :param Kqz: Bulk modulus of grain material in GPa
      :type Kqz: float
      :param Gqz: Shear modulus of grain material in GPa
      :type Gqz: float
      :param Dsh: density the clay
      :type Dsh: float
      :param Ksh: bulk modulus of the clay
      :type Ksh: float
      :param Gsh: shear modulus of the clay
      :type Gsh: float
      :param Dc: density of the cement
      :type Dc: float
      :param Kc: Bulk modulus of cement
      :type Kc: float
      :param Gc: Shear modulus of cement
      :type Gc: float
      :param Db: density of the pore fluid
      :type Db: float
      :param Kb: bulk modulus of the pore fluid
      :type Kb: float
      :param phib_p: adjusted high porosity end memeber
      :type phib_p: float
      :param phic: Critical Porosity
      :type phic: float
      :param sigma: effective stress
      :type sigma: float or array-like
      :param vsh: _description_
      :type vsh: _type_
      :param scheme:
                     Scheme of cement deposition
                             1=cement deposited at grain contacts
                             2=cement deposited at grain surfaces
      :type scheme: int
      :param f:
                reduced shear factor between 0 and 1
                    0=dry pack with inifinitely rough spheres;
                    1=dry pack with infinitely smooth spheres
      :type f: float
      :param Cn: coordination number
      :type Cn: float

      :returns: *array* -- phi,vp1,vp2,vp3,vs1,vs2,vs3: porosity and velocities required for elastic diagnostics bounds















      ..
          !! processed by numpydoc !!

   .. py:method:: screening_plot(phi, vp1, vp2, vp3, cmap=cmap_sand)

      
      plot the rock physics screening crossplot

      :param phi: porosity
      :type phi: array
      :param vp1: lower bound velocities modeled using friable sandstone model by default .
      :type vp1: array
      :param vp2: upper bound velocities modeled using MUSH model by default.
      :type vp2: array
      :param vp3: upper bound velocities modeled using contact cement model by default.
      :type vp3: array
      :param cmap: colormap, can be default colormaps in matplotlib, by default using the customized colormap: cmap_sand
      :type cmap: string, optional

      :returns: *object* -- elastic bounds plot used for rock physics screening















      ..
          !! processed by numpydoc !!

   .. py:method:: normalize(density)
      :staticmethod:

      
      normalize the kde with respect to the maximum value and map value to 0-1
















      ..
          !! processed by numpydoc !!

   .. py:method:: kde_plot(phi, v1, v2, v3, cmap=cmap_kde, vels='Vp', n=64)

      
      plot field data or measurements as 2D probability density functions in the elastic bounds cross plot

      :param phi: porosity
      :type phi: array
      :param vp1: lower bound velocities modeled using friable sandstone model by default .
      :type vp1: array
      :param vp2: upper bound velocities modeled using MUSH model by default.
      :type vp2: array
      :param vp3: upper bound velocities modeled using contact cement model by default.
      :type vp3: array
      :param cmap: colormap, can be default colormaps in matplotlib, by default using the customized colormap: cmap_kde
      :type cmap: string, optional
      :param vels: choose either P wave or S wave velocity for plotting, by default 'Vp'
      :type vels: str, optional
      :param n: grid parameter used in KDE-diffusion, by default 64
      :type n: int, optional

      :returns: *object* -- KDE plot with elastic bounds















      ..
          !! processed by numpydoc !!

   .. py:method:: cst_vels(phi_b, K0, D0, G0, phi, phi_c, Cn, Kc, Gc, Db, Kb, scheme, vsh, Dsh, Dqz, Dc)
      :staticmethod:

      
      compute velocities using constant cement model at different cement amounts

      :param phi_b: adjusted porosity for constant cemnet model
      :type phi_b: float
      :param K0: Bulk modulus of grain material in GPa
      :type K0: float
      :param D0: Density of grain material
      :type D0: float
      :param G0: Shear modulus of grain material in GPa
      :type G0: float
      :param phi: porosity
      :type phi: float
      :param phi_c: Critical Porosity
      :type phi_c: float
      :param Cn: critical porosity
      :type Cn: float
      :param Kc: Bulk modulus of cement
      :type Kc: float
      :param Gc: Shear modulus of cement
      :type Gc: float
      :param Db: density of the pore fluid
      :type Db: float
      :param Kb: bulk modulus of the pore fluid
      :type Kb: float
      :param scheme:
                     Scheme of cement deposition
                             1=cement deposited at grain contacts
                             2=cement deposited at grain surfaces
      :type scheme: int
      :param vsh: shale content
      :type vsh: float
      :param Dsh: density the clay
      :type Dsh: float
      :param Dqz: Density of the grain. not limited to quartz grain
      :type Dqz: float
      :param Dc: density of the cement
      :type Dc: float

      :returns: *array* -- vp,vs: velocities given by constant cement model















      ..
          !! processed by numpydoc !!

   .. py:method:: cst_plot(Dqz, Kqz, Gqz, Dsh, Ksh, Gsh, Dc, Kc, Gc, Db, Kb, phib, phib_p, phi_c, sigma, vsh, Cn, scheme, f)
      :staticmethod:

      
      Diagnostic plot with constant cement model lines

      :param Dqz: Density of the grain. not limited to quartz grain
      :type Dqz: float
      :param Kqz: Bulk modulus of grain material in GPa
      :type Kqz: float
      :param Gqz: Shear modulus of grain material in GPa
      :type Gqz: float
      :param Dsh: density the clay
      :type Dsh: float
      :param Ksh: bulk modulus of the clay
      :type Ksh: float
      :param Gsh: shear modulus of the clay
      :type Gsh: float
      :param Dc: density of the cement
      :type Dc: float
      :param Kc: Bulk modulus of cement
      :type Kc: float
      :param Gc: Shear modulus of cement
      :type Gc: float
      :param Db: density of the pore fluid
      :type Db: float
      :param Kb: bulk modulus of the pore fluid
      :type Kb: float
      :param phib: adjusted high porosity end memeber for constant cement model
      :type phib: float
      :param phib_p: posoities used to drawing the constant cement lines
      :type phib_p: array or list
      :param phi_c: Critical Porosity
      :type phi_c: float
      :param sigma: effective stress
      :type sigma: float
      :param vsh: shale content
      :type vsh: float
      :param Cn: coordination number
      :type Cn: float
      :param scheme:
                     Scheme of cement deposition
                             1=cement deposited at grain contacts
                             2=cement deposited at grain surfaces
      :type scheme: int
      :param f:
                reduced shear factor between 0 and 1
                    0=dry pack with inifinitely rough spheres;
                    1=dry pack with infinitely smooth spheres
      :type f: float

      :returns: *object* -- fig,ax: constant cement diagnostic plot















      ..
          !! processed by numpydoc !!

   .. py:method:: cal_v_const(Dqz, Kqz, Gqz, Dsh, Ksh, Gsh, Dc, Kc, Gc, Db, Kb, phi_b, phi_c, vsh, phi, scheme)
      :staticmethod:

      
      input real data porosity and caculate the theoretical constant cement model velocity value.  Note: input porosity cannot be zero, otherwise the returned velocities are Nan.

      :param Dqz: Density of the grain. not limited to quartz grain
      :type Dqz: float
      :param Kqz: Bulk modulus of grain material in GPa
      :type Kqz: float
      :param Gqz: Shear modulus of grain material in GPa
      :type Gqz: float
      :param Dsh: density the clay
      :type Dsh: float
      :param Ksh: bulk modulus of the clay
      :type Ksh: float
      :param Gsh: shear modulus of the clay
      :type Gsh: float
      :param Dc: density of the cement
      :type Dc: float
      :param Kc: Bulk modulus of cement
      :type Kc: float
      :param Gc: Shear modulus of cement
      :type Gc: float
      :param Db: density of the pore fluid
      :type Db: float
      :param Kb: bulk modulus of the pore fluid
      :type Kb: float
      :param phi_b: _description_
      :type phi_b: _type_
      :param phi_c: Critical Porosity
      :type phi_c: float
      :param vsh: shale content
      :type vsh: float
      :param phi: porosity
      :type phi: array
      :param scheme:
                     Scheme of cement deposition
                             1=cement deposited at grain contacts
                             2=cement deposited at grain surfaces
      :type scheme: int

      :returns: *array* -- vp,vs: constant cement velocities















      ..
          !! processed by numpydoc !!

   .. py:method:: estimate_cem(vcem_seeds, Kqz, Gqz, Ksh, Gsh, phi_c, Cn, Kc, Gc, Db, Kb, scheme, vsh, Dsh, Dqz, Dc)

      
      Estimate cement amount for well log data using constant cement model crossplot.

      :param vcem_seeds: some predefined values at which constant cement lines are calculated
      :type vcem_seeds: array or list
      :param Kqz: Bulk modulus of grain material in GPa
      :type Kqz: float
      :param Gqz: Shear modulus of grain material in GPa
      :type Gqz: float
      :param Ksh: bulk modulus of the clay
      :type Ksh: float
      :param Gsh: shear modulus of the clay
      :type Gsh: float
      :param phi_c: Critical Porosity
      :type phi_c: float
      :param Cn: coordination number
      :type Cn: float
      :param Kc: Bulk modulus of cement
      :type Kc: float
      :param Gc: Shear modulus of cement
      :type Gc: float
      :param Db: density of the pore fluid
      :type Db: float
      :param Kb: bulk modulus of the pore fluid
      :type Kb: float
      :param scheme:
                     Scheme of cement deposition
                             1=cement deposited at grain contacts
                             2=cement deposited at grain surfaces
      :type scheme: int
      :param vsh: shale content
      :type vsh: float
      :param Dsh: density the clay
      :type Dsh: float
      :param Dqz: Density of the grain. not limited to quartz grain
      :type Dqz: float
      :param Dc: density of the cement
      :type Dc: float

      :returns: *array* -- cement amount estimation for each well log data points















      ..
          !! processed by numpydoc !!

   .. py:method:: cement_diag_plot(vcem, Dqz, Kqz, Gqz, Dsh, Ksh, Gsh, Dc, Kc, Gc, Db, Kb, phib, phib_p, phi_c, sigma, vsh, Cn, scheme, f)

      
      _summary_

      :param vcem: _description_
      :type vcem: _type_
      :param Dqz: Density of the grain. not limited to quartz grain
      :type Dqz: float
      :param Kqz: Bulk modulus of grain material in GPa
      :type Kqz: float
      :param Gqz: Shear modulus of grain material in GPa
      :type Gqz: float
      :param Dsh: density the clay
      :type Dsh: float
      :param Ksh: bulk modulus of the clay
      :type Ksh: float
      :param Gsh: shear modulus of the clay
      :type Gsh: float
      :param Dc: density of the cement
      :type Dc: float
      :param Kc: Bulk modulus of cement
      :type Kc: float
      :param Gc: Shear modulus of cement
      :type Gc: float
      :param Db: density of the pore fluid
      :type Db: float
      :param Kb: bulk modulus of the pore fluid
      :type Kb: float
      :param phib: adjusted high porosity end memeber for constant cement model
      :type phib: float
      :param phib_p: posoities used to drawing the constant cement lines
      :type phib_p: array or list
      :param phi_c: Critical Porosity
      :type phi_c: float
      :param sigma: effective stress
      :type sigma: float or array-like
      :param vsh: shale content
      :type vsh: float
      :param Cn: coordination number
      :type Cn: float
      :param scheme:
                     Scheme of cement deposition
                             1=cement deposited at grain contacts
                             2=cement deposited at grain surfaces
      :type scheme: int
      :param f:
                reduced shear factor between 0 and 1
                    0=dry pack with inifinitely rough spheres;
                    1=dry pack with infinitely smooth spheres
      :type f: float

      :returns: *object* -- cross plot for cement estimation















      ..
          !! processed by numpydoc !!

   .. py:method:: plot_rpt(Kdry, Gdry, K0, D0, Kb, Db, Khc, Dhc, phi, sw)
      :staticmethod:

      
      Create RPT plot given computed Impedance and Vp/Vs ratio.

      :param Kdry: effective bulk modulus given by rock physics model
      :type Kdry: float or array
      :param Gdry: effective shear modulus given by rock physics model
      :type Gdry: float or array
      :param K0: bulk modulus of grain
      :type K0: float
      :param D0: density of grain
      :type D0: float
      :param Kb: bulk moduluf of brine
      :type Kb: float
      :param Db: density of brine
      :type Db: float
      :param Khc: bulk modulus of HC
      :type Khc: float
      :param Dhc: density of HC
      :type Dhc: float
      :param phi: porosity
      :type phi: float or array
      :param sw: water saturation
      :type sw: float or array

      :returns: **python onject** (*fig*) -- rpt plot















      ..
          !! processed by numpydoc !!


.. py:class:: utils

   
   Basic calculations for velocities, moduli and stiffness matrix.
















   ..
       !! processed by numpydoc !!
   .. py:method:: V(K, G, rho)
      :staticmethod:

      
      Compute velocity given density and elastic moduli.

      :param K: (GPa): bulk modulus
      :type K: float or array
      :param G: (GPa): shear moulus
      :type G: float or array
      :param rho: (g/m3): density of the frame
      :type rho: float or array

      :returns: *float or array* -- Vp, Vs (m/s): velocity















      ..
          !! processed by numpydoc !!

   .. py:method:: poi(K, G)
      :staticmethod:

      
      Compute poisson's ratio from K an G

      :param K: (GPa): bulk modulus
      :type K: float or array
      :param G: (GPa): shear moulus
      :type G: float or array

      :returns: *float or array* -- Poisson's ratio















      ..
          !! processed by numpydoc !!

   .. py:method:: lame(K, G)
      :staticmethod:

      
      Compute lame constant lamdba from K an G

      :param K: (GPa): bulk modulus
      :type K: float or array
      :param G: (GPa): shear moulus
      :type G: float or array

      :returns: *float or array* -- Poisson's ratio















      ..
          !! processed by numpydoc !!

   .. py:method:: M_from_V(den, vp, vs)
      :staticmethod:

      
      _summary_

      :param den: (g/cm3): bulk density
      :type den: float or array
      :param vp: (m/s): p wave velocity
      :type vp: float or array
      :param vs: (m/s): s wave velocity
      :type vs: float or array

      :returns: *float or array* -- K, G (GPa):bulk and shear moduli















      ..
          !! processed by numpydoc !!

   .. py:method:: write_HTI_matrix(C11, C33, C13, C44, C55)
      :staticmethod:

      
      formulate HTI stiffness matrix

      :param C11: (GPa): stiffness
      :type C11: float
      :param C33: (GPa): stiffness
      :type C33: float
      :param C13: (GPa): stiffness
      :type C13: float
      :param C44: (GPa): stiffness
      :type C44: float
      :param C55: (GPa): stiffness
      :type C55: float

      :returns: *2d array* -- C: 6x6 stiffness matrix















      ..
          !! processed by numpydoc !!

   .. py:method:: write_VTI_compliance(S11, S12, S13, S33, S44)
      :staticmethod:

      
      formulate VTI compliance matrix

      :param S11: (GPa): stiffness
      :type S11: float
      :param S12: (GPa): stiffness
      :type S12: float
      :param S13: (GPa): stiffness
      :type S13: float
      :param S33: (GPa): stiffness
      :type S33: float
      :param S44: (GPa): stiffness
      :type S44: float

      :returns: *2d array* -- S: 6x6 compliance matrix















      ..
          !! processed by numpydoc !!

   .. py:method:: write_VTI_matrix(C11, C33, C13, C44, C66)
      :staticmethod:

      
      formulate VTI stiffness matrix

      :param C11: (GPa): stiffness
      :type C11: float
      :param C33: (GPa): stiffness
      :type C33: float
      :param C13: (GPa): stiffness
      :type C13: float
      :param C44: (GPa): stiffness
      :type C44: float
      :param C66: (GPa): stiffness
      :type C66: float

      :returns: *2d array* -- C: 6x6 stiffness matrix















      ..
          !! processed by numpydoc !!

   .. py:method:: write_matrix(C11, C22, C33, C12, C13, C23, C44, C55, C66)
      :staticmethod:

      
      formulate general 6x6 stiffness matrix in Voigt notation

      :param C11: (GPa): stiffness
      :type C11: float
      :param C22: (GPa): stiffness
      :type C22: float
      :param C33: (GPa): stiffness
      :type C33: float
      :param C12: (GPa): stiffness
      :type C12: float
      :param C13: (GPa): stiffness
      :type C13: float
      :param C23: (GPa): stiffness
      :type C23: float
      :param C44: (GPa): stiffness
      :type C44: float
      :param C55: (GPa): stiffness
      :type C55: float
      :param C66: (GPa): stiffness
      :type C66: float

      :returns: *2d array* -- C: 6x6 stiffness matrix















      ..
          !! processed by numpydoc !!

   .. py:method:: write_iso(K, G)
      :staticmethod:

      
      formulate isotropic 6x6 stiffness matrix in Voigt notation

      :param K:
      :type K: float or array

          (GPa): bulk modulus
      G : float or array
          (GPa): shear moulus

      :returns: *2d array* -- C: 6x6 stiffness matrix















      ..
          !! processed by numpydoc !!

   .. py:method:: crack_por(crd, alpha)
      :staticmethod:

      
      compute crack porosity from crack aspect ratio and crack density

      :param crd: (unitless): crack density
      :type crd: float or array
      :param alpha: crack aspect ratio
      :type alpha: float or array

      :returns: *float or array* -- cpor (frac): crack porosity















      ..
          !! processed by numpydoc !!

   .. py:method:: v_to_c_VTI(Vp0, Vp45, Vp90, Vs0, Vsh90, den)
      :staticmethod:

      
      compute stiffness matrix given velocity measurements along different directions

      :param Vp0: (km/s): incident angle dependent velocity measurements
      :type Vp0: float or array
      :param Vp45: (km/s): incident angle dependent velocity measurements
      :type Vp45: float or array
      :param Vp90: (km/s): incident angle dependent velocity measurements
      :type Vp90: float or array
      :param Vs0: (km/s): incident angle dependent velocity measurements
      :type Vs0: float or array
      :param Vsh90: (km/s): incident angle dependent velocity measurements
      :type Vsh90: float or array
      :param den: (g/cm3):density of the sample
      :type den: float or array

      :returns: *2d array* -- C: VTI stiffness matrix















      ..
          !! processed by numpydoc !!


