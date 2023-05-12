:py:mod:`rockphypy.AVO`
=======================

.. py:module:: rockphypy.AVO


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.AVO.AVO




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

   .. py:method:: Aki_Richard(theta, vp1, vp2, vs1, vs2, den1, den2)
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


