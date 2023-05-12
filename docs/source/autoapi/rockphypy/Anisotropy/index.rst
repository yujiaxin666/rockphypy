:py:mod:`rockphypy.Anisotropy`
==============================

.. py:module:: rockphypy.Anisotropy


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.Anisotropy.Anisotropy




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


