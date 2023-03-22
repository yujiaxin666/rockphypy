"""
Critical porosity model and Gassmann Fluid substitution
=======================================================

This pages gives a simple examples of using the classes ``EM`` and ``Fluid`` of ``rockphypy`` to compute Critical porosity model and perform Gassmann Fluid substitution. 
"""
# %% 

# Importing aux libraries
import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams['font.size']=14
plt.rcParams['font.family']='arial'

# %% 
import rockphypy # import the module for rock physics
from rockphypy import EM # import the "effective medium" EM module 
# import the 'Fluid' module 
from rockphypy import Fluid


# %% 
#
# Critical porosity model
# ~~~~~~~~~~~~~~~~~~~~~~~
# 
# Remember Nur’s hypothesis: There is a critical (structure-dependent)
# porosity at which the framework stiffness goes to zero. The simple yet 
# powerful critical porosity model is defined as: 
# 
# .. math::
#         K_{\text {dry }}=K_{0}\left(1-\frac{\phi}{\phi_{\mathrm{c}}}\right)
#
# .. math::
#         \mu_{\text {dry }}=\mu_{0}\left(1-\frac{\phi}{\phi_{\mathrm{c}}}\right)
# 
# where :math:`K_0` and :math:`\mu_0` are the mineral bulk and shear moduli.
# 

# %% 
#
# Gassmann Fluid substitution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute dry frame moduli using rock physics models such as critical porosity 
# model is straightforward as shown above. To compute saturated elastic moduli, Gassmann's explicit equations for fluid substition can be applied. Gassmann's relations not only allows us to perform fluid 
# substition (when one fluid is replaced
# with another), but to predict saturated-rock moduli from dry-rock moduli, and # vice versa. 
# 
# Here we show how to compute saturated moduli via Gassmann thoery. 
# 
# The equation we are gonna use is 
# 
# .. math::
#         K_{\text {sat }}=K_{\text {dry }}+\frac{\left(1-K_{\text {dry }} / K_{0}\right)^{2}}{\phi / K_{\mathrm{fl}}+(1-\phi) / K_{0}-K_{\text {dry }} / K_{0}^{2}}
# 
# 
# This is one of the equivelent expressions of the well known GS-Biot theory for fluid substitution:
# 
# .. math::
#         \frac{K_{\text {sat }}}{K_{0}-K_{\text {sat }}}=\frac{K_{\text {dry }}}{K_{0}-K_{\text {dry }}}+\frac{K_{\text {fl }}}{\phi\left(K_{0}-K_{\text {fl }}\right)}, \quad \mu_{\text {sat }}=\mu_{\text {dry }}
# 
#
# Now Let's firstly estimate the effective dry-rock and water saturated moduli # using critical porosity model. Assume that the rock is Arenite with 40% 
# critical porosity.
#

# %% 

# specify model parameters
phic=0.4
phi=np.linspace(0.001,phic,100,endpoint=True) # solid volume fraction = 1-phi
K0, G0= 37,44
Kw = 2.2
Kg = 0.5

# %% 

# Compute dry-rock moduli
K_dry, G_dry= EM.cripor(K0, G0, phi, phic)
# saturate rock with water 
Ksat, Gsat = Fluid.Gassmann(K_dry,G_dry,K0,Kw,phi)

# %% 

# plot
# sphinx_gallery_thumbnail_number = 1
plt.figure(figsize=(6,6))
plt.xlabel('Porosity')
plt.ylabel('Bulk modulus [GPa]')
plt.title('V, R, VRH, HS bounds')
plt.plot(phi, K_dry,label='dry rock K')
plt.plot(phi, Ksat,label='saturated K')

plt.legend(loc='best')
plt.grid(ls='--')

# %% 
# 
# We can see from the figure that effective bulk modulus increases when the 
# rock is saturated.
# 
