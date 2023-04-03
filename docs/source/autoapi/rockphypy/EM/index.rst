:py:mod:`rockphypy.EM`
======================

.. py:module:: rockphypy.EM


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.EM.EM




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

      :returns: *_type_* -- C_eff: effective moduli of cracked, transversely isotropic rocks

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

      :returns: *_type_* -- C_eff: effective moduli in 6x6 matrix form.















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
      :type Ks: float
      :param Gs: Shear modulus of matrix in GPa
      :type Gs: float
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
      :type crd: float

      :returns: *float* -- K_dry,G_dry: dry elastic moduli of cracked medium















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
      :param f: _descripvolume fraction of inclusion phase tion_
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


