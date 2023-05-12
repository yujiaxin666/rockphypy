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
# The question is as follows: Calculate the self-consistent effective bulk and shear moduli, :math:`K_{eff}^{SC}` and :math:`G_{eff}^{SC}`, for a water-saturated rock consisting of spherical quartz grains (aspect ratio :math:`\alpha` = 1) and total porosity 0.3. The pore space consists of spherical pores :math:`\alpha` = 1 and thin penny-shaped cracks (:math:`\alpha` = 1e^−2). The thin cracks have a porosity of 0.01, whereas the remaining porosity (0.29) is made up of the spherical pores.
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

#%%
# The effective moduli of the two phase composite as a function of the volume fraction of the soft fluid filled crack can be calculated as follow:  
#

# two phase K=[K1,K2]
K=[37,2.25]
# G=[G1,G2]
G=[44,0,]
Alpha=[1,0.1]

frac = np.linspace(0,1,50)
K_eff = np.zeros(frac.size)
G_eff = np.zeros(frac.size)
for i, val in enumerate(frac):
    X=[1-val, val]
    K_eff[i],G_eff[i]= EM.Berryman_sc(K,G,X,Alpha)

#%%

# sphinx figure 
# sphinx_gallery_thumbnail_number = 1
plt.figure(figsize=(5,5))
plt.plot(frac,K_eff,'-k',lw=3,label='K_eff')
plt.plot(frac,G_eff,'-b',lw=2,label='G_eff')
plt.xlabel('Volume fraction of soft phase')
plt.ylabel('Effective modulus')
plt.legend()
# %%
# As can be shown in the figure, the effective shear modulus of the two phase composite becomes 0 when the volume fraction of the soft phase is approximately 45% for an aspect ratio of 0.1. this prediction is very similar to the critical porosity model which predicts a suspension of grain in the fluid when the porosity exceeds about 0.4. 
# 
# However, feel free to change the aspect ratio for the soft phase from 0.1 to 0.01, then the modelling results of the effective shear modulus becomes zero when the volume fraction of the soft phase is approximately 12%. 
#