"""
Simple anisotropy
=================
This exampl shows how to use the ``rockphypy.Anisotropy`` class to compute Thomsen parameters, Weak anistropic approximation of phase velocities and Backus averaging.
"""

# %%


import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams['font.size']=14
plt.rcParams['font.family']='arial'


# %%

import rockphypy # import the module
from rockphypy import Anisotropy


# %%
#
# Notice that stressed-induced seismic anisotropy: if an initially isotropic
# material develops anisotropy due to applied stress, the anisotropy must have at least orthorhombic symmetry. 
# 
# Thomsen parameters
# ~~~~~~~~~~~~~~~~~~
# The Elastic constant  of a transversely isotropic elastic material  in terms of Voigt notation can be represented as 
# 
# .. math::
#         \begin{bmatrix}
#           c_{11} & c_{12} & c_{13} & 0 & 0 & 0 \\
#           c_{12} & c_{11} & c_{13} & 0 & 0 & 0 \\
#           c_{13} & c_{13} & c_{33} & 0 & 0 & 0 \\
#            0 & 0 & 0 & c_{44} & 0 & 0 \\
#            0 & 0 & 0 & 0 & c_{44} & 0\\
#            0 & 0 & 0 & 0 & 0 & c_{66}
#          \end{bmatrix}
# 
# where 
#
# .. math::
#       c_{66}=\frac{1}{2}\left(c_{11}-c_{12}\right)
# 
# 
# Here z-axis is the unique symmetry axis, isotropic in xy-plane.
# 
# Thomsen (1986) simplified the elasticity of a VTI material by introducing the anisotropy parameters
# 
# P wave anisotropy: 
#
# .. math::
#       \varepsilon=\frac{c_{11}-c_{33}}{2 c_{33}}
#  
# S wave anisotropy:
#
# .. math::
#       \gamma=\frac{c_{66}-c_{44}}{2 c_{44}}
#   
# 
# Moveout parameter: 
#
# .. math::
#       \delta=\frac{\left(c_{13}+c_{44}\right)^{2}-\left(c_{33}-c_{44}\right)^{2}}{2 c_{33}\left(c_{33}-c_{44}\right)}
#   
# 
# P wave velocity: 
#
# .. math::
#       \alpha=\sqrt{c_{33} / \rho}
#  
# 
# S wave velocity: 
#
# .. math::
#       \beta=\sqrt{c_{44} / \rho}
# 
# 
# Intepretation of :math:`\varepsilon`, :math:`\gamma` and :math:`\delta`: 
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The preceding equations are valid for any strength of VTI anisotropy, since they are just definitions. 
# 
# But for weak anisotropy, :math:`\varepsilon` is usually called "P-wave anisotropy", as it can be seen to approximately describe the fractional difference between the P-wave velocities parallel and orthogonal to the symmetry axis
# 
# .. math::
#       \varepsilon \approx \frac{V_{\mathrm{P}}\left(90^{\circ}\right)-V_{\mathrm{P}}\left(0^{\circ}\right)}{V_{\mathrm{P}}\left(0^{\circ}\right)}
# 
# 
# the constant :math:`\gamma` can be seen to describe the fractional difference between the SH-wave velocities parallel and orthogonal to the symmetry axis, which is equivalent to the difference between the velocities of S-waves polarized parallel and normal to the symmetry axis, both propagating normal to the symmetry axis:
# 
# .. math::
#       \gamma \approx \frac{V_{\mathrm{SH}}\left(90^{\mathrm{o}}\right)-V_{\mathrm{SV}}\left(90^{\circ}\right)}{V_{\mathrm{SV}}\left(90^{\circ}\right)}=\frac{V_{\mathrm{SH}}\left(90^{\mathrm{o}}\right)-V_{\mathrm{SH}}\left(0^{\mathrm{o}}\right)}{V_{\mathrm{SH}}\left(0^{\circ}\right)}
# 
# 
# :math:`\delta` is called moveout parameter as the small-offset normal moveout (NMO) velocity is affected by VTI anisotropy, :math:`\delta` goes into the equation for NMO velocities, :math:`V_{NMO,P}`, :math:`V_{NMO,SV}`, and :math:`V_{NMO,SH}` for
# P-, SV-, and SH-modes calculation: 
# 
# .. math::
#       V_{\mathrm{NMO}, \mathrm{P}}=\alpha \sqrt{1+2 \delta}
# 
# 
# .. math::
#       V_{\mathrm{NMO}, \mathrm{SV}}=\beta \sqrt{1+2 \sigma}, \quad \sigma=\left(\frac{\alpha}{\beta}\right)^{2}(\varepsilon-\delta)
# 
# 
# .. math::
#       V_{\mathrm{NMO}, \mathrm{SH}}=\beta \sqrt{1+2 \gamma}
# 
# 
# In terms of the Thomsen parameters, the three phase velocities for weak anisotropy can
# be approximated as
# 
# .. math::
#       V_{\mathrm{P}}(\theta) \approx \alpha\left(1+\delta \sin ^{2} \theta \cos ^{2} \theta+\varepsilon \sin ^{4} \theta\right)
# 
# 
# .. math::
#       V_{\mathrm{SV}}(\theta) \approx \beta\left[1+\frac{\alpha^{2}}{\beta^{2}}(\varepsilon-\delta) \sin ^{2} \theta \cos ^{2} \theta\right]
# 
# 
# .. math::
#       V_{\mathrm{SH}}(\theta) \approx \beta\left(1+\gamma \sin ^{2} \theta\right)
# 
#
# Backus average
# ~~~~~~~~~~~~~~
# The Backus average is used to model a finely stratified medium as a single homogeneous medium.all materials are linearly elastic; there are no sources of intrinsic energy dissipation, such as friction or viscosity;
# and the layer thickness must be much smaller than the seismic wavelength. 
# 
# For a periodically layered medium with isotropic layers of :math:`n` materials {:math:`n_1`, :math:`n_2`, ...} having concentrations {:math:`f_1`,  :math:`f_2`,....}, :math:`f_1` + :math:`f_2` +...=1 and elastic moduli (Lamé coefficients) :math:`\lambda_1`, :math:`G_1`, :math:`\lambda_2`, :math:`G_2`,... :math:`\lambda_n`, :math:`G_n`, the effective anisotropy of the layered medium is given by
# 
# .. math::
#       C_{11}=\left\langle\frac{4 \mu(\lambda+\mu)}{\lambda+2 \mu}\right\rangle+\left\langle\frac{1}{\lambda+2 \mu}\right\rangle^{-1}\left\langle\frac{\lambda}{\lambda+2 \mu}\right\rangle^{2}
# 
# 
# .. math::
#       C_{33}=\left\langle\frac{1}{\lambda+2 \mu}\right\rangle^{-1}
# 
# 
# .. math::
#       C_{13}=\left\langle\frac{1}{\lambda+2 \mu}\right\rangle^{-1}\left\langle\frac{\lambda}{\lambda+2 \mu}\right\rangle
# 
# 
# .. math::
#       C_{44}=\left\langle\frac{1}{\mu}\right\rangle^{-1}
# 
# 
# .. math::
#       C_{66}=\langle\mu\rangle
# 
# 
# .. math::
#       C_{12}=C_{11}-2C_{66}
# 
# 
#
# Examples
# ~~~~~~~~
# Let's estimate angle dependent weak anisotropy of a layered medium using backus average model.
#


# %%


# specify model parameters
lamda1, G1=5,5
lamda2, G2=1,1
den1=2.25
den2=2.0
V1,V2= 0.5,0.5 # volumetric fraction
# Compute anisotropic elastic moduli
V = [V1,V2]
lamda = [lamda1,lamda2]
G=[G1,G2]
C11,C33,C13,C44,C66= Anisotropy.Backus(V,lamda,G)
# compute effective density
den= np.dot(V,np.array([den1,den2]))
# compute angle dependent anisotropy from layering
theta=np.linspace(0,90,50)
VP, VSV, VSH,epsilon,gamma,delta = Anisotropy.Thomsen(C11,C33,C13,C44,C66,den, theta)


# %%


# plot
plt.figure(figsize=(6,6))
plt.xlabel('angle')
plt.ylabel('velocity')
plt.title('Anisotropy from layering')
plt.plot(theta, VP, label='VqP')
plt.plot(theta, VSV,label='VSH')
plt.plot(theta, VSH,label='VqSV')
plt.legend(loc='best')


# %%
# **Reference** : Mavko, G., Mukerji, T. and Dvorkin, J., 2020. The rock physics handbook. Cambridge university press.
# 
