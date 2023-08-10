.. rockphypy documentation master file, created by
   sphinx-quickstart on Fri Mar  3 16:42:58 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

About
=====
An extensive Python library for rock physics modelling
******************************************************


Installation
------------

``rockphypy`` is available through `PYPI <https://pypi.org/project/rockphypy/>`_ , and may be installed using ``pip``: 


   ``$ pip install rockphypy``



Overview
--------
``rockphypy`` is a Python-based **open-source library** dedicated to rock physics modeling. It caters to individuals interested in rock physics and empowers them to apply rock physics tools in their work and research effectively.

Built upon the foundations of `Stanford SRB tools  <https://github.com/StanfordRockPhysics/SRBToolbox>`_ , rockphypy has successfully migrated and optimized a wide range of rock physics models from Stanford SRB Matlab tools to Python. Additionally, it extends the functionality by introducing new functions and practical workflows. Notably, ``rockphypy`` offers a diverse set of valuable modules such as ``Anisotropy``, ``AVO``, ``BW`` (Batzle&Wang), ``EM`` (Effective Medium), ``Empirical``, ``Fluid``, ``GM`` (Granular Medium), ``Permeability``, ``QI`` (Quantitative Intepretation) and ``utils``. For a more detailed understanding of each module, please refer to the API Documentation.

From the ground up, ``rockphypy`` is designed to be easily embedded in various workflows i.e, machine learning, probabilistic frameworks, bayesian inversion that leverage the rich resources of python open-source ecosystem.


Contribution and Development
----------------------------
See the timeline of ``rockphypy`` `here <https://user-images.githubusercontent.com/45630390/259761842-5689968c-7683-41e4-864a-0dca791a38a0.png>`_  for more details about the development.

Driven by a commitment to constant growth and development, rockphypy welcomes contributions from the community. Whether through `pull request <https://github.com/yujiaxin666/RockPhysicsPy/pulls>`_  or by raising questions and making requests via `issue creation <https://github.com/yujiaxin666/RockPhysicsPy/issues>`_. Your active involvement is highly encouraged. Together, we can foster a collaborative and thriving environment for rock physics enthusiasts.



Example code
------------
Below is a simple example demonstrating the usage of the :class:`~rockphypy.EM` and :class:`~rockphypy.Fluid` from the ``rockphypy`` library to create saturated elastic bounds. For further examples, please explore the **Rock physics basics** and **Applications** tabs on the left.


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

The motivation behind the creation of this Python library stems from the Isabella's personal experience and passion for rock physics. As an avid fan of the Python language and the open-source community, Isabella embarked on her PhD journey in rock physics.

Isabella encountered challenges in finding a user-friendly rock physics library to assist with their head-drain exercise when she began her PhD study at NTNU. This experience sparked the realization of a need for an accessible and convenient rock physics tool.

Simultaneously, Isabella's fascination with the theoretical versatility of rock physics modeling grew, leading she to delve into classical papers and books on the subject. The rock physics handbook and other excellent books and papers (see **literature**) played a significant role in shaping her understanding of rock physics.

Fortuitously, Isabella had the opportunity to visit Stanford during her PhD, where Isabella combined her passion for rock physics with her love for Python and open-source development. This led to the creation of this Python library, serving as the final piece of their PhD work.

The library aims to provide an intuitive and practical solution for rock physics enthusiasts, drawing from the author's firsthand experience and knowledge gained through PhD research.


Table of contents
==================
.. currentmodule:: rockphypy
.. toctree::
   
   getting_started/index
   advanced_examples/index
   autoapi/index
   errata
   
Contents:

.. toctree::
   :maxdepth: 3

Indices
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

