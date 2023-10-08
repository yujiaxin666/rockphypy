"""
Varying Patchiness Cement Model(VPCM)
=====================================
"""
# %%

import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams['font.size']=14
plt.rcParams['font.family']='arial'

# %%

#import rockphypy # import the module 
from rockphypy import GM, EM



# %%
# The Hertz-Mindlin thoery indicates that unconsolidated sands exhibit stress sensitivity that arises from grain-grain contact and friction-resisted rotation. The widely applied contact cement model (Dvorkin and Nur, 1996) implies that when the sands are cemented, the resulting sandstone will become completely stress insensitive. However, in real cases, cemented sandstone can still possess stress sensitivity that might originate from the inhomogeneous spatial distribution of the cement within the grain packing. Futhermore, laboratory tests have shown that the sandstone shows significant stress sensitivity when subjected to stress release.  Yu et al. (2023) present a new rock physics model to quantitatively describe the stress sensitivity of weakly to moderately cemented sandstone during effective stress release. The model is built upon the patchy cement model (PCM) and incorporates microscopic observations of cement cracking and crumbling. To account for the reduced coherence of the cement coherence upon stress release, a cement diluting factor :math:`\alpha`, which helps analyze stress sensitivity changes during stress removal is introduced. 
#
# Below we performed a thorough analysis of the VPCM. 
#
# To understand the idea of varying patchiness cement model, we can start by looking at the stress sensitivity curves of disconnected patchy cement model and connected patchy cement model computed at using the same porosity and same amount of cement. Then the only difference between the two computed sandstone models (connected patchy cement sandstone and disconnected patch cement sandstone) is the microgeometry assumed. Let's model the stress sensitivity curves for porosity=0.36 using PCM:

# %%

# specify model parameters
Dqz, Kqz, Gqz = 2.65, 36, 42 ## grain density, bulk and shear modulus 
Dsh, Ksh, Gsh = 2.7, 21, 7 # shale/clay density, bulk and shear modulus
Dc,Kc, Gc =2.65, 36, 42 # cement density, bulk and shear modulus
vsh=0 # clay fraction 

_,_,K0=EM.VRH(np.array([vsh,1-vsh]),np.array([Ksh,Kqz])) # clay fraction can be considered.
_,_,G0=EM.VRH(np.array([vsh,1-vsh]),np.array([Gsh,Gqz]))

phic=0.4 # critical porosity

phi=np.linspace(1e-6,phic,100)
Cn=6 # coordination number
v_cem=0.1 # critical cement limit 
v_ci=0.111
scheme=2 # cement deposition 
f_=0.5 #reduce shear factor
f=0.8 # effective cement fraction assumed in PCM
phi=0.36 # porosity 

sigma=np.linspace(1e-7,20,100)

# connected patchy cement 
Kdry1,Gdry1=GM.pcm(f,sigma, K0,G0,phi, phic,v_cem,v_ci, Kc,Gc,Cn, 'stiff',scheme,f_)
# disconnected patchy cement 
Kdry2,Gdry2=GM.pcm(f,sigma, K0,G0,phi, phic,v_cem,v_ci, Kc,Gc,Cn, 'soft',scheme,f_)

#plot
fig=plt.figure(figsize=(6,6))
plt.xlabel('Peff MPa')
plt.ylabel('$K_{dry}$ GPa')
plt.xlim(0,20)
plt.ylim(2.5,5.5)
plt.title('PCM at $\phi$=0.36')
plt.plot(sigma,Kdry1,'-r',lw=3,label='Connected patchy')
plt.plot(sigma,Kdry2,'-b',lw=3,label='Disconnected patchy')
plt.fill_between(sigma, Kdry1, Kdry2, color='grey', alpha=0.5)
plt.legend(loc='best')

# %%
# Between these two curves, there is a noticeable gap, and VPCM seeks to clarify what it represents within that range.
#


# %%
#
# Build Varying Patchiness Cement Model 
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Recall that the patchy cement model represents effective medium comprising a binary mixture of cemented sandstone and unconsolidated sand. At high porosity end, the two end member elasticitities :math:`K_{unc}` and :math:`K_{cem}` are modelled with and Walton theory with slip factor and contact cement model respectively. 
# 
# .. math::
#         K_{connected}=K_{cem}+\frac{(1-fcc)}{\left(K_{unc}-K_{cem}\right)^{-1}+fcc\left(K_{cem}+\frac{4}{3} \mu_{cem}\right)^{-1}}
# 
# .. math::
#         K_{disconnected}=K_{unc}+\frac{fcc}{\left(K_{cem}-K_{unc}\right)^{-1}+(1-fcc)\left(K_{unc}+\frac{4}{3} \mu_{unc}\right)^{-1}}
# 
# Microcracks can develop in the grain contact cement when sandstone is subjected to effective stress release. How to introduce the crack is the big challenge. VPCM address this challenge by cement diluting: the connected patchy sandstone is being replaced by disconnected patchy cement sandstone upon stress release. Then the elasticities given by VPCM at high-porosity end member are:
#
# .. math::
#           K_{vpcm} = K_{connected}(1-\alpha)+ \alpha K_{disconnected} 
#
# The diluting factor :math:`\alpha` quantifies how much of connected patchy cement has been replaced by disconnected patchy cement upon stress unloading. It is allowed to be stressdependent by
#
# .. math::
#          \alpha=\left(1-\frac{\sigma^{\prime}}{\sigma_0^{\prime}}\right)^m
#
# The effective dry rock moduli at smaller porosity are computed using Lower Hashin-Strikmann mixing as done in Soft Sand model and patchy cement model. 
# 

# %%
# The VPCM has been implemented in the ``GM`` module by calling the GM.vpcm. for the forward modeling, predefined the diluting schedule using ``GM.diluting``. 
#

# %%
k=1 # full diluting no cement crushing 
sigma0 = sigma.max()
m= 1 
alpha = GM.diluting(k,sigma0,sigma,m)

Kdry4,Gdry4= GM.vpcm(alpha, f,sigma,K0,G0,phi, phic, v_cem,v_ci, Kc,Gc, Cn,scheme,f_)


# sphinx figure 
# sphinx_gallery_thumbnail_number = 2

fig=plt.figure(figsize=(6,6))
plt.xlabel('Effective stress (MPa)')
plt.ylabel(r'$K_{\rm dry}$ (GPa)')
plt.xlim(0,20)
plt.ylim(2.5,5.5)
#plt.title('(b) VPCM at $\phi$=0.36, f=0.8')
plt.plot(sigma,Kdry1,'-r',lw=3,label='Connected patchy cement model')
plt.plot(sigma,Kdry2,'-b',lw=3,label='Disconnected patchy cement model')
#plt.fill_between(sigma, Kdry1, Kdry2, color='grey', alpha=0.5)
plt.plot(sigma,Kdry4,'--k',lw=3,label='VPCM')
plt.legend(loc='best')
#plt.text(12,4.8,'$\\lambda$:[0,1]')



# %%
# **Reference**: 
#
# - Yu, J., Duffaut, K., & Avseth, P. (2023). Stress sensitivity of elastic moduli in high-porosity cemented sandstone—Heuristic models and experimental data. Geophysics, 88(4), MR185-MR194.
#  
