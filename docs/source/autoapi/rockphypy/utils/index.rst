:py:mod:`rockphypy.utils`
=========================

.. py:module:: rockphypy.utils


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.utils.utils




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


