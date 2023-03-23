:py:mod:`rockphypy.AVO`
=======================

.. py:module:: rockphypy.AVO

.. autoapi-nested-parse::

   @File    :   AVO.py
   @Time    :   2023/02/12 16:36:09
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

   rockphypy.AVO.AVO




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


