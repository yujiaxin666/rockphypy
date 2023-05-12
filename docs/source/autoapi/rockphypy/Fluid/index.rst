:py:mod:`rockphypy.Fluid`
=========================

.. py:module:: rockphypy.Fluid


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.Fluid.Fluid




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


