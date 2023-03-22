"""
Berryman Self consistent approach
=================================
Berryman Self consistent approach for calculating the elastic moduli of N phase composite
"""

# %%

import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams['font.size']=14
plt.rcParams['font.family']='arial'


# %%

import rockphypy # import the module 
from rockphypy import EM 
from rockphypy import Fluid

# %%
# Berryman's self-consistent approximations for N-phase composites are described by the equations: 
# 
#.. math::
#       \sum_{\mathrm{i}=1}^N x_{\mathrm{i}}\left(K_{\mathrm{i}}-K_{\mathrm{SC}}^{e f f}\right) P^{e f f \mathrm{i}}=0 
#
#.. math::
#        \sum_{\mathrm{i}=1}^N x_{\mathrm{i}}\left(\mu_{\mathrm{i}}-\mu_{\mathrm{SC}}^{e f f}\right) Q^{e f f \mathrm{i}}=0
# 
#

# %%
# The question is as follows: Calculate the self-consistent effective bulk and shear moduli, :math:`K_{eff}^{SC}`and :math:`G_{eff}^{SC}`, for a water-saturated rock consisting of spherical quartz grains (aspect ratio :math:`\alpha` = 1) and total porosity 0.3. The pore space consists of spherical pores :math:`\alpha` = 1 and thin penny-shaped cracks (:math:`\alpha` = 1e^−2). The thin cracks have a porosity of 0.01, whereas the remaining porosity (0.29) is made up of the spherical pores.
# 
#
# There are three phases in the composite, i.e. quartz grain, water filled spherical pore, and water filled thin cracks. We can use the method ``EM.Berryman_sc`` to easily solve this exercise to find the effective moduli of the composite.
#


# %%

# three phase K=[K1,K2,K3]
K=[37,2.25,2.25]
# G=[G1,G2,G3]
G=[44,0,0]
# X=[frac1, frac2, frac3]
X=[0.7,0.29,0.01]
# Alpha= [alpha1, alpha2, alpha3]
Alpha=[1,1,0.01]

# effective moduli
K_sc,G_sc= EM.Berryman_sc(K,G,X,Alpha)


# %%


print('K_eff and G_eff of the composite are {:.2f} GPa and {:.2f} GPa, respectively'.format(K_sc,G_sc))

