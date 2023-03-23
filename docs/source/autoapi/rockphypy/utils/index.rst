:py:mod:`rockphypy.utils`
=========================

.. py:module:: rockphypy.utils

.. autoapi-nested-parse::

   @File    :   utils.py
   @Time    :   2023/01/16 12:06:24
   @Author  :   Jiaxin Yu
   @Contact :   yujiaxin666@outlook.com
   @License :   (C)Copyright 2020-2021, Jiaxin Yu

   ..
       !! processed by numpydoc !!


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


