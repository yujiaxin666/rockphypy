.. rockphypy documentation master file, created by
   sphinx-quickstart on Fri Mar  3 16:42:58 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to rockphypy's documentation!

About
=====
Open-source rock physic toolboxs in Python.
*******************************************

Overview
--------
``rockphypy`` is a Python-based **open-source rock physics modelling library** created for people who are interested in Rock physics and want to applied rock physics tools in their work and research.

The library refactor most of rock physics models in `Stanford SRB tools  <https://github.com/StanfordRockPhysics/SRBToolbox>`_ from Matlab to Python. 
In addition, more functions and various practical workflows are created. ``rockphypy`` provides a bunch of usefel classes, i.e. `Anisotropy`, `AVO`, `BW (Batzle&Wang)`,`EM(Effective medium)`,`Empirical`, `Fluid`,`GM(Granular Medium)`,`Permeability`and `utils`. See the API Documentation for further details. 

From the ground up, ``rockphypy`` is designed to be easily embedded in various workflows i.e, machine learning, probabilistic frameworks, bayesian inversion that leverage the rich resources of python open-source. 

Since a constant state of development and growth is aimed, please feel free to contribute by `pulling request <https://github.com/yujiaxin666/RockPhysicsPy/pulls>`_ or pose your questions and requests by `creating issues <https://github.com/yujiaxin666/RockPhysicsPy/issues>`_. 



Installation
------------

``rockphypy`` is available through ``PyPI`` and may be installed using ``pip``: ::


   $ pip install rockphypy

Example code
------------
Here's an example showing the usage of :class:`~rockphypy.EM.EM` and :class:`~rockphypy.Fluid.Fluid` for creating saturated elastic bounds. 

.. plot::
   :include-source:

   from rockphypy import EM
   from rockphypy import Fluid

   # specify model parameters
   phi=np.linspace(0,1,100,endpoint=True) # solid volume fraction = 1-phi
   K0, G0= 37,44 # moduli of grain material
   Kw, Gw= 2.2,0 # moduli of water 
   # VRH bounds
   volumes= np.vstack((1-phi,phi)).T
   M= np.array([K0,Kw])
   K_v,K_r,K_h=EM.VRH(volumes,M)
   # Hashin-Strikmann bound 
   K_UHS,G_UHS= EM.HS(1-phi, K0, Kw,G0,Gw, bound='upper')
   # Critical porosity model
   phic=0.4 # Critical porosity
   phi_=np.linspace(0.001,phic,100,endpoint=True) # solid volume fraction = 1-phi
   K_dry, G_dry= EM.cripor(K0, G0, phi_, phic)# Compute dry-rock moduli
   Ksat, Gsat = Fluid.Gassmann(K_dry,G_dry,K0,Kw,phi_)# saturate rock with water

   # plot
   plt.figure(figsize=(6,6))
   plt.xlabel('Porosity')
   plt.ylabel('Bulk modulus [GPa]')
   plt.title('V, R, VRH, HS bounds')
   plt.plot(phi, K_v,label='K Voigt')
   plt.plot(phi, K_r,label='K Reuss = K HS-')
   plt.plot(phi, K_h,label='K VRH')
   plt.plot(phi, K_UHS,label='K HS+')
   plt.plot(phi_, Ksat,label='K CriPor')
   plt.legend(loc='best')
   plt.grid(ls='--')


Motivation 
----------

The motivation comes from author's personal experience in rock physics study and the fact that she is a big fans of python and open source community. 
When she began her phd study working on rock physics, she's looking for an easy to use rock physics library for finishing her head drain exercise for reservoir seismic lectures. 
On the other hand, the author is facinated by the theoretical versaility of rock physics modelling and submerge herself into reading classical books for rock physics and pratical applications.  The rock physics handbook indeed gave here inspriations and help her build a solid understanding of the rock physics. eventually the author did her visit at stanford and finishing the python library as the last piece of her phd work. this is indeed a beautiful experience. 



Table of contents
==================
.. currentmodule:: rockphypy
.. toctree::
   
   getting_started/index
   autoapi/index

Contents:

.. toctree::
   :maxdepth: 3

Indices
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

