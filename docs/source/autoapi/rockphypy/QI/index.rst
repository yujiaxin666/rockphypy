:py:mod:`rockphypy.QI`
======================

.. py:module:: rockphypy.QI


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.QI.QI




Attributes
~~~~~~~~~~

.. autoapisummary::

   rockphypy.QI.cmap_sand
   rockphypy.QI.cmap_kde
   rockphypy.QI.cmap_cem


.. py:data:: cmap_sand

   

.. py:data:: cmap_kde

   

.. py:data:: cmap_cem

   

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


