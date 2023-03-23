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
   utils/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.Anisotropy
   rockphypy.AVO
   rockphypy.BW
   rockphypy.utils
   rockphypy.GM
   rockphypy.EM




.. py:class:: Anisotropy

   
   Effective models, coordinate transform, anisotropic parameters and phase velocities that can be applied to anisotropic media.
















   ..
       !! processed by numpydoc !!
   .. py:method:: Thomsen(C33, C13, C44, C66, den, theta)

      
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

   .. py:method:: Thomsen_Tsvankin(C22, C33, C12, C13, C23, C44, C55, C66)

      
      Elastic constants of an orthorhombic elastic medium defined by Tsvankin’s notation for weak elastic anisotropy assuming the vertical symmetry axis is along the x3 direction.

      :param Cij: Stiffnesses in GPa
      :type Cij: float

      :returns: *floats* -- Thomsen-Tsvankin parameters















      ..
          !! processed by numpydoc !!

   .. py:method:: Backus(lamda, G)

      
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

   .. py:method:: Backus_log(Vs, Den, Depth)

      
      Computes Backus Average from log data, notice that the Depth is 1d Vector including each top depth of layer and also the bottom of last layer.

      :param Vp: P wave velocities of layers [Vp1,Vp2...Vpn], size N
      :type Vp: array
      :param Vs: S wave velocities of layers [Vs1,Vs2...Vsn], size N
      :type Vs: array
      :param Den: Densities of layers, size N
      :type Den: array
      :param Depth: 1d depth, ATTENTION: each depth point corresponds to the top of thin isotropic layer, the bottom of the sedimentary package is the last depth point. [dep1,dep2,,,,,,depn, depn+1], size N+1
      :type Depth: array

      :returns: *array-like* -- Stiffness coeffs and averaged density















      ..
          !! processed by numpydoc !!

   .. py:method:: vel_azi_HTI(Den, azimuth)

      
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

   .. py:method:: vel_azi_VTI(Den, azimuth)

      
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

   .. py:method:: Bond_trans(theta, axis=3)

      
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
   .. py:method:: AVO_HTI(D2, C1, C2, theta, azimuth)

      
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
      :type theta: float or array-like
      :param azimuth: azimuth angle, in degree unit
      :type azimuth: float or array-like

      :returns: *float or array-like* -- PP reflectivities















      ..
          !! processed by numpydoc !!

   .. py:method:: Aki_Richard(vp1, vp2, vs1, vs2, den1, den2)

      
      Aki-Richard approximation to PP reflectivity.

      :param theta: incident angle, degree
      :type theta: float or array-like
      :param vp1: P wave velocity of layer 1, m/s
      :type vp1: float or array-like
      :param vp2: P wave velocity of layer 2, m/s
      :type vp2: float or array-like
      :param vs1: S wave velocity of layer 1, m/s
      :type vs1: float or array-like
      :param vs2: S wave velocity of layer 2, m/s
      :type vs2: float or array-like
      :param den1: density of layer 1, kg/m3
      :type den1: float or array-like
      :param den2: density of layer 2, kg/m3
      :type den2: float or array-like

      :returns: *float or array-like* -- R_pp: P wave reflectivity
                R_ps: PS reflectivity
                Rpp0: intercept
                gradient















      ..
          !! processed by numpydoc !!

   .. py:method:: zoeppritz(vs1, rho1, vp2, vs2, rho2, theta)

      
      Reflection & Transmission coefficients calculated using full Zoeppritz equations.

      :param vp1: P wave velocity of layer 1, m/s
      :type vp1: float or array-like
      :param vs1: S wave velocity of layer 1, m/s
      :type vs1: float or array-like
      :param rho1: density of layer 1, kg/m3
      :type rho1: float or array-like
      :param vp2: P wave velocity of layer 2, m/s
      :type vp2: float or array-like
      :param vs2: S wave velocity of layer 2, m/s
      :type vs2: float or array-like
      :param rho2: density of layer 2, kg/m3
      :type rho2: float or array-like
      :param theta: incident angle, degree
      :type theta: float or array-like

      :returns: *float or array-like* -- Rpp,Rps: PP and PS reflectivity















      ..
          !! processed by numpydoc !!

   .. py:method:: AVO_abe(vs1, d1, vp2, vs2, d2)

      
      Copied from RPT matlab tools func: avo_abe

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

   .. py:method:: EI_ref(Vs, rho, theta, SP, norm=True)

      
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
      :param theta: incident angle
      :type theta: float or array-like
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

   .. py:method:: AVO_ortho(b1, e11, d11, e12, d12, g1, rho1, a2, b2, e21, d21, e22, d22, g2, rho2, the)

      
      calculates the reflectivity in the symmetry plane for interfaces between 2 orthorhombic media

      :param a1: _description_
      :type a1: _type_
      :param b1: _description_
      :type b1: _type_
      :param e11: _description_
      :type e11: _type_
      :param d11: _description_
      :type d11: _type_
      :param e12: _description_
      :type e12: _type_
      :param d12: _description_
      :type d12: _type_
      :param g1: _description_
      :type g1: _type_
      :param rho1: _description_
      :type rho1: _type_
      :param a2: _description_
      :type a2: _type_
      :param b2: _description_
      :type b2: _type_
      :param e21: _description_
      :type e21: _type_
      :param d21: _description_
      :type d21: _type_
      :param e22: _description_
      :type e22: _type_
      :param d22: _description_
      :type d22: _type_
      :param g2: _description_
      :type g2: _type_
      :param rho2: _description_
      :type rho2: _type_
      :param the: _description_
      :type the: _type_

      :returns: *_type_* -- _description_















      ..
          !! processed by numpydoc !!


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


.. py:class:: utils

   
   Basic calculations for velocities, moduli and stiffness matrix.
















   ..
       !! processed by numpydoc !!
   .. py:method:: V(G, rho)

      
      Compute velocity given density and elastic moduli.
      :param K: bulk modulus
      :type K: GPa
      :param G: shear moulus
      :type G: GPa
      :param rho: density of the frame
      :type rho: g/m3

      :returns: *Vp, Vs (m/s)* -- velocity

      Written by Jiaxin Yu (July 2021)















      ..
          !! processed by numpydoc !!

   .. py:method:: poi(G)

      
      Compute poisson's ratio from K an G

      :param K: bulk modulus
      :type K: GPa
      :param G: shear modulus
      :type G: GPa

      :returns: *nu* -- Poisson's ratio

      Written by Jiaxin Yu (July 2021)















      ..
          !! processed by numpydoc !!

   .. py:method:: lame(G)

      
      Compute lame constant lamdba from K an G

      :param K: bulk modulus
      :type K: GPa
      :param G: shear modulus
      :type G: GPa

      :returns: *nu* -- Poisson's ratio

      Written by Jiaxin Yu (July 2021)















      ..
          !! processed by numpydoc !!

   .. py:method:: M_from_V(vp, vs)

      
      Compute K and G from velocities and density

      :param den: bulk density
      :type den: g/cm3
      :param vp: p wave velocity
      :type vp: m/s
      :param vs: s wave velocity
      :type vs: m/s

      :returns: *K, G (GPa)* -- bulk and shear moduli

      Written by Jiaxin Yu (July 2021)















      ..
          !! processed by numpydoc !!

   .. py:method:: write_HTI_matrix(C33, C13, C44, C55)

      
      formulate HTI stiffness matrix

      :param C11: stiffness
      :type C11: GPa
      :param C13: stiffness
      :type C13: GPa
      :param C23: stiffness
      :type C23: GPa
      :param C33: stiffness
      :type C33: GPa
      :param C44: stiffness
      :type C44: GPa
      :param C55: stiffness
      :type C55: GPa

      :returns: *C* -- 6x6 stiffness matrix















      ..
          !! processed by numpydoc !!

   .. py:method:: write_VTI_compliance(S12, S13, S33, S44)

      
      formulate VTI compliance matrix

      :param S11: compliance
      :type S11: GPa
      :param S12: compliance
      :type S12: GPa
      :param S13: compliance
      :type S13: GPa
      :param S33: compliance
      :type S33: GPa
      :param S44: compliance
      :type S44: GPa

      :returns: *_type_* -- _description_















      ..
          !! processed by numpydoc !!

   .. py:method:: write_VTI_matrix(C33, C13, C44, C66)

      
      formulate VTI stiffness matrix

      :param C11: stiffness
      :type C11: GPa
      :param C33: stiffness
      :type C33: GPa
      :param C13: stiffness
      :type C13: GPa
      :param C44: stiffness
      :type C44: GPa
      :param C65: stiffness
      :type C65: GPa

      :returns: *C* -- 6x6 stiffness matrix















      ..
          !! processed by numpydoc !!

   .. py:method:: write_matrix(C22, C33, C12, C13, C23, C44, C55, C66)

      
      formulate general 6x6 stiffness matrix in Voigt notation

      :param Cij: stiffness
      :type Cij: GPa

      :returns: *C* -- 6x6 stiffness matrix















      ..
          !! processed by numpydoc !!

   .. py:method:: write_iso(G)

      
      formulate isotropic 6x6 stiffness matrix in Voigt notation

      :param Cij: stiffness
      :type Cij: GPa

      :returns: *C* -- 6x6 stiffness matrix















      ..
          !! processed by numpydoc !!

   .. py:method:: crack_por(alpha)

      
      compute crack porosity from crack aspect ratio and crack density

      :param crd: crack density
      :type crd: unitless
      :param alpha: crack aspect ratio
      :type alpha: unitless

      :returns: *cpor (frac)* -- crack porosity















      ..
          !! processed by numpydoc !!

   .. py:method:: v_to_c_VTI(Vp45, Vp90, Vs0, Vsh90, den)

      
      _summary_

      :param Vp0: indident angle dependent velocity measurements
      :type Vp0: km/s
      :param Vp45: indident angle dependent velocity measurements
      :type Vp45: km/s
      :param Vp90: indident angle dependent velocity measurements
      :type Vp90: km/s
      :param Vs0: indident angle dependent velocity measurements
      :type Vs0: km/s
      :param Vsh90: indident angle dependent velocity measurements
      :type Vsh90: km/s
      :param den: density of the sample
      :type den: g/cm3

      :returns: *C* -- VTI stiffness matrix















      ..
          !! processed by numpydoc !!


.. py:class:: GM

   
   Contact based granular medium models and extensions.
















   ..
       !! processed by numpydoc !!
   .. py:method:: ThomasStieber(phi_sh, vsh)

      
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

   .. py:method:: silty_shale(Kq, Gq, Ksh, Gsh)

      
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

   .. py:method:: shaly_sand(C, Kss, Gss, Kcc, Gcc)

      
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

   .. py:method:: contactcement(G0, Kc, Gc, phi, phic, Cn, scheme)

      
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

   .. py:method:: hertzmindlin(G0, phic, Cn, sigma, f)

      
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

   .. py:method:: softsand(G0, phi, phic, Cn, sigma, f)

      
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

   .. py:method:: Walton(G0, phic, Cn, sigma, f)

      
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

   .. py:method:: johnson(G0, n, phi, epsilon, epsilon_axial, path='together')

      
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

   .. py:method:: stiffsand(G0, phi, phic, Cn, sigma, f)

      
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

   .. py:method:: constantcement(K0, G0, Kc, Gc, phi, phic, Cn, scheme)

      
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

   .. py:method:: MUHS(G0, Kc, Gc, phi, phi_b, phic, Cn, scheme)

      
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

   .. py:method:: Digby(G0, phi, Cn, sigma, a_R)

      
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

   .. py:method:: pcm(sigma, K0, G0, phi, phic, v_cem, v_ci, Kc, Gc, Cn, mode, scheme, f_)

      
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


.. py:class:: EM

   
   classical bounds and inclusion models.
















   ..
       !! processed by numpydoc !!
   .. py:method:: VRH(M)

      
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

   .. py:method:: cripor(G0, phi, phic)

      
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

   .. py:method:: cripor_reuss(Mf, phic, den=False)

      
      In the suspension domain, the effective bulk and shear moduli of the rock can be estimated by using the Reuss (isostress) average.

      :param M0: The solid phase modulus or density
      :type M0: float or array-like
      :param Mf: The pore filled phase modulus or density
      :type Mf: float or array-like
      :param phic: critical porosity
      :type phic: float
      :param den: If False: compute the reuss average for effective modulus of two mixing phases. If true, compute avearge density using mass balance, which corresponds to voigt average. Defaults to False.
      :type den: bool, optional

      :returns: *float or array-like* -- M (GPa/g.cc): average modulus or average density

      .. rubric:: References

      - Section 7.1 Rock physics handbook 2nd edition















      ..
          !! processed by numpydoc !!

   .. py:method:: HS(K1, K2, G1, G2, bound='upper')

      
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

   .. py:method:: Eshelby_Cheng(G, phi, alpha, Kf, mat=False)

      
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

      :returns: *_type_* -- C_eff: effective moduli of cracked, transversely isotropic rocks

      .. rubric:: References

      - section 4.14 in The Rock Physics Handbook















      ..
          !! processed by numpydoc !!

   .. py:method:: hudson(G, Ki, Gi, alpha, crd, order=1, axis=3)

      
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

      :returns: *_type_* -- C_eff: effective moduli in 6x6 matrix form.















      ..
          !! processed by numpydoc !!

   .. py:method:: hudson_rand(G, Ki, Gi, alpha, crd)

      
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

   .. py:method:: hudson_ortho(G, Ki, Gi, alpha, crd)

      
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

   .. py:method:: hudson_cone(G, Ki, Gi, alpha, crd, theta)

      
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

   .. py:method:: Berryman_sc(G, X, Alpha)

      
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

   .. py:method:: PQ_vectorize(Gm, Ki, Gi, alpha)

      
      compute geometric strain concentration factors P and Q for prolate and oblate spheroids according to Berymann (1980).See also: Berryman_sc, Berryman_func

      :param Km: Shear modulus of matrix phase. For Berryman SC       approach, this corresponds to the effective moduli of the composite.
      :type Km: float
      :param Gm: Bulk modulus of matrix phase. For Berryman SC approach, this corresponds to the effective moduli of the composite.
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

   .. py:method:: Berryman_func(K, G, X, Alpha)

      
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

   .. py:method:: Swiss_cheese(Gs, phi)

      
      Compute effective elastic moduli via "Swiss cheese" model with spherical pores. "Swiss cheese" model assumes a dilute distribution of spherical inclusions embedded in an * *unbounded* * homogenous solid.  It takes the "noninteracting assumption" in which all cavities (pores) are independent so that their contributions can be added.

      :param Ks: Bulk modulus of matrix in GPa
      :type Ks: float
      :param Gs: Shear modulus of matrix in GPa
      :type Gs: float
      :param phi: porosity
      :type phi: float or array-like

      :returns: *float or array-like* -- Kdry,Gdry (GPa): effective elastic moduli















      ..
          !! processed by numpydoc !!

   .. py:method:: SC(Ks, Gs, iter_n)

      
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

   .. py:method:: Dilute_crack(Gs, cd)

      
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

   .. py:method:: OConnell_Budiansky(G0, crd)

      
      O’Connell and Budiansky (1974) presented equations for effective bulk and shear moduli of a cracked medium with randomly oriented dry penny-shaped cracks (in the limiting case when the aspect ratio α goes to 0)

      :param K0: bulk modulus of background medium
      :type K0: float
      :param G0: shear modulus of background medium
      :type G0: float
      :param crd: crack density
      :type crd: float

      :returns: *float* -- K_dry,G_dry: dry elastic moduli of cracked medium















      ..
          !! processed by numpydoc !!

   .. py:method:: OConnell_Budiansky_fl(G0, Kfl, crd, alpha)

      
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

   .. py:method:: OC_R_funcs(crd, nu_0, w)

      
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

   .. py:method:: PQ(Gm, Ki, Gi, alpha)

      
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

   .. py:method:: DEM(t, params)

      
      ODE solver tutorial: https://physics.nyu.edu/pine/pymanual/html/chap9/chap9_scipy.html.
















      ..
          !! processed by numpydoc !!

   .. py:method:: Berryman_DEM(Gm, Ki, Gi, alpha, phi)

      
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

