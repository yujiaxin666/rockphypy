:py:mod:`rockphypy.QI`
======================

.. py:module:: rockphypy.QI


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.QI.QI




.. py:class:: QI

   
   Useful functionalities for quantitative intepretation and relevant tasks.
















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


