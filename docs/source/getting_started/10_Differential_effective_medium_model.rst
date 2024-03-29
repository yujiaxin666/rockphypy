
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "getting_started\10_Differential_effective_medium_model.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        :ref:`Go to the end <sphx_glr_download_getting_started_10_Differential_effective_medium_model.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_getting_started_10_Differential_effective_medium_model.py:


Differential Effective Medium (DEM)
===================================

.. GENERATED FROM PYTHON SOURCE LINES 8-15

.. code-block:: python3


    import numpy as np 
    import matplotlib.pyplot as plt
    plt.rcParams['font.size']=14
    plt.rcParams['font.family']='arial'









.. GENERATED FROM PYTHON SOURCE LINES 16-21

.. code-block:: python3


    import rockphypy # import the module 
    from rockphypy import EM 
    from rockphypy import Fluid








.. GENERATED FROM PYTHON SOURCE LINES 22-123

Berryman's Differential effective medium (DEM) model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The DEM theory models two-phase composites by incrementally adding inclusions of one phase (phase 2) to the matrix phase. The matrix begins as phase 1 (when the concentration of phase 2 is zero) and is changed at each step as a new increment of phase 2 material is added. The process is continued until the desired proportion of the constituents is reached. The process of incrementally adding inclusions to the matrix is really a thought experiment and should not be taken to provide an accurate description of the true evolution of rock porosity in nature.

The coupled system of ordinary differential equations for the effective bulk and
shear moduli, :math:`K^{DEM}_{eff}` and :math:`\mu^{DEM}_{eff}`, respectively, are (Berryman, 1992b)

.. math::
      (1-y) \frac{d}{d y}\left[K_{D E M}^{e f f}(y)\right]=\left(K_{2}-K_{D E M}^{e f f}\right) P^{e f f 2}(y)


.. math::
      (1-y) \frac{d}{d y}\left[\mu_{D E M}^{e f f}(y)\right]=\left(\mu_{2}-\mu_{D E M}^{e f f}\right) Q^{e f f 2}(y)


with initial conditions  :math:`K^{DEM}_{eff} (0)=K_1` and :math:`\mu^{DEM}_{eff}(0)=\mu_1`, where :math:`K_1` and :math:`\mu_1` are the bulk and shear moduli of the initial host material (phase 1), :math:`K_2` and :math:`\mu_2` are the bulk and shear moduli of the incrementally added inclusions (phase 2), and y is the concentration of phase 2.
For fluid inclusions and voids, y equals the porosity.

Expressions for the volumetric and deviatoric strain concentration factors are:

.. math::
      P^{m n}=\frac{1}{3}T_{i i j j}^{m n}

.. math::
      Q^{m n}=\frac{1}{5}(T_{i j i j}^{m n}-T_{i i j j}^{m n})

.. math::
      T_{i i j j}^{m n}=\frac{3 F_{1}}{F_{2}}

.. math::
      T_{i j i j}^{m n}-\frac{1}{3} T_{i i j j}^{m n}=\frac{2}{F_{3}}+\frac{1}{F_{4}}+\frac{F_{4} F_{5}+F_{6} F_{7}-F_{8} F_{9}}{F_{2} F_{4}}

where

.. math::
      F_{1}=1+A\left[\frac{3}{2}(f+\theta)-R\left(\frac{3}{2} f+\frac{5}{2} \theta-\frac{4}{3}\right)\right]

.. math::
      F_{2}=1+A\left[1+\frac{3}{2}(f+\theta)-\frac{1}{2} R(3 f+5 \theta)\right]+B(3-4 R)+\frac{1}{2} A(A+3 B)(3-4 R)\left[f+\theta-R\left(f-\theta+2 \theta^{2}\right)\right]  \\

.. math::
      F_{3}=1+A\left[1-\left(f+\frac{3}{2} \theta\right)+R(f+\theta)\right]

.. math::
      F_{4}=1+\frac{1}{4} A[f+3 \theta-R(f-\theta)] 

.. math::
      F_{5}=A\left[-f+R\left(f+\theta-\frac{4}{3}\right)\right]+B \theta(3-4 R) 

.. math::
      F_{6}=1+A[1+f-R(f+\theta)]+B(1-\theta)(3-4 R) 

.. math::
      F_{7}=2+\frac{1}{4} A[3 f+9 \theta-R(3 f+5 \theta)]+B \theta(3-4 R) 

.. math::
      F_{8}=A\left[1-2 R+\frac{1}{2} f(R-1)+\frac{1}{2} \theta(5 R-3)\right]+B(1-\theta)(3-4 R) 

.. math::
      F_{9}=A[(R-1) f-R \theta]+B \theta(3-4 R)

with A, B, and R given by

.. math::
      A=\mu_{\mathrm{i}} / \mu_{\mathrm{m}}-1

.. math::
      B=\frac{1}{3}\left(\frac{K_{\mathrm{i}}}{K_{\mathrm{m}}}-\frac{\mu_{\mathrm{i}}}{\mu_{\mathrm{m}}}\right)

.. math::
      R=\frac{\left(1-2 v_{\mathrm{m}}\right)}{2\left(1-v_{\mathrm{m}}\right)}

The functions :math:`\theta` and :math:`f` are given by 

.. math::
      \theta=\begin{cases}
      \{\frac{\alpha}{\left(\alpha^{2}-1\right)^{3 / 2}}\left[\alpha\left(\alpha^{2}-1\right)^{1 / 2}-\cosh ^{-1} \alpha\right]\\
      \frac{\alpha}{\left(1-\alpha^{2}\right)^{3 / 2}}\left[\cos ^{-1}\alpha-\alpha\left(1-\alpha^{2}\right)^{1 / 2}\right]
      \end{cases}

for prolate and oblate spheroids, respectively, and

.. math::
      f=\frac{\alpha^{2}}{1-\alpha^{2}}(3 \theta-2)


Note that :math:`\alpha <1` for oblate spheroids and :math:`\alpha >1`  for prolate spheroids

For spherical pores:

.. math::
      P=\frac{K_{\mathrm{m}}+\frac{4}{3} \mu_{\mathrm{m}}}{K_{\mathrm{i}}+\frac{4}{3} \mu_{\mathrm{m}}} 


.. math::
      Q=\frac{\mu_{\mathrm{m}}+\zeta_{\mathrm{m}}}{\mu_{\mathrm{i}}+\zeta_{\mathrm{m}}}

Fluid Effect:
^^^^^^^^^^^^^ 
Dry cavities can be modeled by setting the inclusion moduli to zero. Fluid-saturated cavities are simulated by setting the inclusion shear modulus to zero. 
%%

.. GENERATED FROM PYTHON SOURCE LINES 123-139

.. code-block:: python3



    # DEM modelling the crack 
    Gi =0
    Ki =0 # dry pore
    alpha = 1
    # Initial values
    K_eff0 = 37   # host mineral bulk modulus 
    G_eff0 = 45   # host mineral shear modulus 

    # Make time array for solution, desired fraction of inclusion 
    tStop = 1

    K_dry_dem, G_dry_dem,t= EM.Berryman_DEM(K_eff0,G_eff0, Ki, Gi, alpha,tStop)









.. GENERATED FROM PYTHON SOURCE LINES 140-151

.. code-block:: python3



    # plot
    plt.figure(figsize=(6,6))
    plt.xlabel('Porosity')
    plt.ylabel('Bulk modulus [GPa]')
    plt.title('DEM dry-pore modelling')
    plt.plot(t, K_dry_dem)
    plt.grid(ls='--')





.. image-sg:: /getting_started/images/sphx_glr_10_Differential_effective_medium_model_001.png
   :alt: DEM dry-pore modelling
   :srcset: /getting_started/images/sphx_glr_10_Differential_effective_medium_model_001.png
   :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 152-154

**Reference**: Mavko, G., Mukerji, T. and Dvorkin, J., 2020. The rock physics handbook. Cambridge university press.



.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  0.087 seconds)


.. _sphx_glr_download_getting_started_10_Differential_effective_medium_model.py:

.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-example




    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download Python source code: 10_Differential_effective_medium_model.py <10_Differential_effective_medium_model.py>`

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download Jupyter notebook: 10_Differential_effective_medium_model.ipynb <10_Differential_effective_medium_model.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
