"""
Digby model
===========
"""

# %%

import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams['font.size']=14
plt.rcParams['font.family']='arial'
plt.rcParams['axes.labelpad'] = 10.0
plt.rcParams["figure.figsize"] = (6,6)


# %%


import rockphypy # import the module 
from rockphypy import GM
from rockphypy import  utils

# %%
# The Digby model gives effective moduli for a dry, random packing of identical elastic spherical particles. Neighboring particles are initially firmly bonded across small, flat, circular regions of radius a. Outside these adhesion surfaces, the shape of each particle is assumed to be ideally smooth (with a continuous first derivative). Notice that this condition differs from that of Hertz, where the shape of a particle is not smooth at the intersection of the spherical surface and the plane of contact. Digby’s normal and shear stiffnesses under hydrostatic pressure :math:`P` are (Digby, 1981)
#
# .. math::
#       S_{n}=\frac{4 \mu b}{1-v}, \quad S_{\tau}=\frac{8 \mu a}{2-v}
# 
# Parameter b can be found from the relation
# 
# .. math::
#       \frac{b}{R}=\left[d^{2}+\left(\frac{a}{R}\right)^{2}\right]^{1 / 2}
# 
# where :math:`d` satisfies the cubic equation
# 
# .. math::
#       d^{3}+\frac{3}{2}\left(\frac{a}{R}\right)^{2} d-\frac{3 \pi(1-v) P}{2 C(1-\phi) \mu}=0
# 
# Like other effective granular medium models, the Digby's model is only strictly valid for seismic waves of sufficiently low frequency those for which
# 
# .. math::
#       \rho \omega^{2} R^{2} /(\lambda+2 \mu)<\rho \omega^{2} R^{2} / \mu \ll 1.
# 
# Example
# ^^^^^^^
#

# %%


a_R=np.array([0., 0.02, 0.05, 0.07,0.5])
G0= 44
K0=37
D0= 2.65
phi= 0.36
D= (1-phi)* D0
Cn=9
sigma=np.linspace(1e-6,20,20)# Mpa


K_dry=np.zeros((a_R.size, sigma.size)) # each row corresponds to a fixed a/R ratio
G_dry=np.zeros((a_R.size, sigma.size)) # each row corresponds to a fixed a/R ratio

# notice that the np.roots can only accept rank 1 array, so when sigma is a array with multiple entries, a for loop is used.
for i in range(a_R.size):
    for j, val in enumerate(sigma):

        Keff, Geff= GM.Digby(K0, G0, phi, Cn, val, a_R[i] )

        K_dry[i,j]= Keff
        G_dry[i,j]= Geff

Vp, Vs= utils.V(K_dry, G_dry, D)


# %%


plt.figure(figsize=(8,7))
plt.plot(sigma, Vp[0,:],'ko-', label='a/R=0',clip_on=False)
plt.plot(sigma, Vp[1,:],'k+-',label='a/R=0.02',clip_on=False)
plt.plot(sigma, Vp[2,:],'kD-',label='a/R=0.05',clip_on=False)

plt.plot(sigma, Vs[0,:],'bo-',clip_on=False)
plt.plot(sigma, Vs[1,:],'b+-',clip_on=False)
plt.plot(sigma, Vs[2,:],'bD-',clip_on=False)
plt.xlim(0,20)
plt.ylim(0,2000)
plt.xlabel('Pressure')
plt.ylabel('Velocities')
plt.legend()


# %%

plt.figure(figsize=(8,7))
plt.plot(sigma, Vp[0,:]/Vs[0,:],'ko-', label='a/R=0',clip_on=False)
plt.plot(sigma, Vp[1,:]/Vs[1,:],'k+-',label='a/R=0.02',clip_on=False)
plt.plot(sigma, Vp[2,:]/Vs[2,:],'kD-',label='a/R=0.05',clip_on=False)
plt.plot(sigma, Vp[3,:]/Vs[3,:],'k^-', label='a/R=0.07',clip_on=False)
plt.plot(sigma, Vp[4,:]/Vs[4,:],'k*-',label='a/R=0.5',clip_on=False)


plt.xlim(0,20)
#plt.ylim(0,2)
plt.xlabel('Pressure')
plt.ylabel('Vp/Vs')
plt.legend(loc='best', bbox_to_anchor=(0.06,0.5, 1., .102)) #  (x, y, width, height) 

# %%
# This indicates that when the bounding radius decreases, the Vp/Vs ratio will increase. when the ratio of boundin raidus to the grain radius increases, which can be used to describe the increasing contact cement saturation, we can see the stress dependency of the Vp/Vs ratio will vanishes. 
#  
