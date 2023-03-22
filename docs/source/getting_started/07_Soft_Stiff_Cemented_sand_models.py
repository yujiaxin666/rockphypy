"""
Contact based models
====================
"""

# %%

import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams['font.size']=14
plt.rcParams['font.family']='arial'


# %%

import rockphypy # import the module 
from rockphypy import GM


# %%
#
# Soft-Sand model
# ~~~~~~~~~~~~~~~
# In Walton and Hertz–Mindlin contact models, the porosity is fixed at a critical porosity :math:`\phi_c \approx 0.4` for randomly densed packed spheres of equal size. To find the effective moduli of unconsolidated sand at a different porosity :math:`\phi`, a heuristic modified Hashin–Shtrikman lower bound is used as an interpolator between the
# Hertz–Mindlin moduli :math:`K_{\mathrm{HM}}`  and :math:`\mu_{\mathrm{HM}}` at porosity :math:`\phi_c` and the solid grain moduli :math:`K` and :math:`\mu` at zero porosity, the elastic moduli at different porosity can be computed as:
# 
# .. math::
#       K_{\mathrm{eff}}=\left[\frac{\phi / \phi_{c}}{K_{\mathrm{HM}}+\frac{4}{3} \mu_{\mathrm{HM}}}+\frac{1-\phi / \phi_{c}}{K+\frac{4}{3} \mu_{\mathrm{HM}}}\right]^{-1}-\frac{4}{3} \mu_{\mathrm{HM}}
# 
# 
# .. math::
#       \mu_{\mathrm{eff}}= \left[\frac{\phi / \phi_{c}}{\mu_{\mathrm{HM}}+\frac{\mu_{\mathrm{HM}}}{6}\left(\frac{9 K_{\mathrm{HM}}+8 \mu_{\mathrm{HM}}}{K_{\mathrm{HM}}+2 \mu_{\mathrm{HM}}}\right)}+\frac{1-\phi / \phi_{c}}{\mu+\frac{\mu_{\mathrm{HM}}}{6}\left(\frac{9 K_{\mathrm{HM}}+8 \mu_{\mathrm{HM}}}{K_{\mathrm{HM}}+2 \mu_{\mathrm{HM}}}\right)}\right]^{-1} -\frac{\mu_{\mathrm{HM}}}{6}\left(\frac{9 K_{\mathrm{HM}}+8 \mu_{\mathrm{HM}}}{K_{\mathrm{HM}}+2 \mu_{\mathrm{HM}}}\right) 
# 
# 
# where :math:`K` and :math:`\mu` are bulk and shear moduli of grain material, repectively. 
# 
# Soft-sand model is also called unconsolidated sand model or friable sand model e.g. in Avseth et al. (2010).
#
# Stiff-Sand model
# ~~~~~~~~~~~~~~~~
# A counterpart to the soft-sand model is the “stiff-sand” model, which uses precisely the same end-members as in the soft-sand model but connects them with a heuristic modified Hashin–Shtrikman upper bound as a stiff interpolator: 
# 
# .. math::
#       K_{\mathrm{eff}}=\left[\frac{\phi / \phi_{c}}{K_{\mathrm{HM}}+\frac{4}{3} \mu}+\frac{1-\phi / \phi_{c}}{K+\frac{4}{3} \mu}\right]^{-1}-\frac{4}{3} \mu
# 
# 
# .. math::
#       \mu_{\mathrm{eff}}=   \left[\frac{\phi / \phi_{c}}{\mu_{\mathrm{HM}}+\frac{\mu}{6}\left(\frac{9 K+8 \mu}{K+2 \mu}\right)}+\frac{1-\phi / \phi_{c}}{\mu+\frac{\mu}{6}\left(\frac{9 K+8 \mu}{K+2 \mu}\right)}\right]^{-1}-\frac{\mu}{6}\left(\frac{9 K+8 \mu}{K+2 \mu}\right) 
# 
#
# Dvorkin’s Cemented-Sand Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Soft- and Stiff-Sand models describe how the elastic moduli of unconsolidated sand vary in the porosity-stress plane. Dvorkin’s Cemented-Sand Model instead allows us to compute the moduli of sand pack in which the cement deposits at the *grain contact*. As a result, this model is also called **Contact cement model**. The contact cement fills the crack-like spaces near the grain contacts. This has the effect of very rapidly stiffening the rock
# with very little change in porosity. This cement tends to eliminate further sensitivity. In this model, the effective bulk modulus :math:`K_{dry}` and shear modulus :math:`G_{dry}` of dry rock are:
# 
# .. math::
#       K_{\mathrm{eff}}=\frac{1}{6} C\left(1-\phi_{0}\right) M_{\mathrm{c}} \hat{S}_{\mathrm{n}}
# 
# .. math::
#       \mu_{\mathrm{eff}}=\frac{3}{5} K_{\mathrm{eff}}+\frac{3}{20} C\left(1-\phi_{0}\right) \mu_{\mathrm{c}} \hat{S}_{\tau} 
# 
# .. math::
#       M_c= K_c+ \frac{4\mu_c}{3}
# 
# where :math:`\phi_c` is critical porosity; :math:`K_s` and :math:`\mu_s` are the bulk and shear moduli of the grain material, respectively; :math:`K_c` and  :math:`\mu_c` are the bulk and shear moduli of the cement material, respectively; :math:`M_c` is the compressional modulus of the cement; and n is the coordination number. The parameters :math:`\hat{S}_{\mathrm{n}} `and   :math:`\hat{S}_{\tau}` are:
# 
# .. math::
#       \hat{S}_{\mathrm{n}}=A_{\mathrm{n}} \alpha^{2}+B_{\mathrm{n}} \alpha+C_{\mathrm{n}}
#
# .. math::
#       A_{\mathrm{n}}=-0.024153 \Lambda_{\mathrm{n}}^{-1.3646}
#
# .. math::
#       B_{\mathrm{n}}=0.20405 \Lambda_{\mathrm{n}}^{-0.89008}
#
# .. math::
#       C_{\mathrm{n}}=0.00024649 \Lambda_{\mathrm{n}}^{-1.9864}
#
# .. math::
#       \hat{S}_{\tau}=A_{\tau} \alpha^{2}+B_{\tau} \alpha+C_{\tau}
#
# .. math::
#       A_{\tau}=-10^{-2}\left(2.26 v^{2}+2.07 v+2.3\right) \Lambda_{\tau}^{0.079 v^{2}+0.1754 v-1.342}
#
# .. math::
#       B_{\tau}=\left(0.0573 v^{2}+0.0937 v+0.202\right) \Lambda_{\tau}^{0.0274 v^{2}+0.0529 v-0.8765}
#
# .. math::
#       C_{\tau}=10^{-4}\left(9.654 v^{2}+4.945 v+3.1\right) \Lambda_{\tau}^{0.01867 v^{2}+0.4011 v-1.8186}
#
# .. math::
#       \Lambda_{\mathrm{n}}=\frac{2 \mu_{\mathrm{c}}}{\pi \mu} \frac{(1-v)\left(1-v_{\mathrm{c}}\right)}{\left(1-2 v_{\mathrm{c}}\right)}
#
# .. math::
#       \Lambda_{\tau}=\frac{\mu_{\mathrm{c}}}{\pi \mu}
#
# .. math::
#       \alpha=\frac{a}{R}
# 
# By assuming that porosity reduction in sands is due to cementation
# only, we can relate the parameter \alpha to the current porosity of cemented sand \phi. For Scheme 1 in which all cement is deposited at grain contacts
# 
# .. math::
#       \alpha=2\left[\frac{\phi_{0}-\phi}{3 C\left(1-\phi_{0}\right)}\right]^{1 / 4}=2\left[\frac{S \phi_{0}}{3 C\left(1-\phi_{0}\right)}\right]^{1 / 4}
# 
# 
# For scheme 2, in which cement is evenly deposited on the grain surface:
# 
# .. math::
#       \alpha=\left[\frac{2\left(\phi_{0}-\phi\right)}{3\left(1-\phi_{0}\right)}\right]^{1 / 2}=\left[\frac{2 S \phi_{0}}{3\left(1-\phi_{0}\right)}\right]^{1 / 2}
# 
# 
# 
# From The handbook of rock physics (Mavko, 2020)
#
# Constant cement model 
# ~~~~~~~~~~~~~~~~~~~~~
# As introduced by Avseth et al. (2000), Constant cement model assumes that sands of varying porosity all have the same amount of contact cement. Porosity variation is solely due to non-contact pore-filling material (e.g., deteriorating sorting). This model is contact cement model blend with soft-sand model at an adjusted high porosity end memeber :math:`\phi_b`. Firstly, porosity reduces from the initial sand-pack porosity to porosity
# :math:`\phi_b`, dry-rock bulk and shear moduli at this porosity (:math:`K_b` and :math:`\mu_b`, respectively) are calculated
# from the contact-cement model. Then the dry-rock bulk :math:`K_{dry}` and shear :math:`\mu_{dry}`
# moduli at a smaller porosity :math:`\phi` are then interpolated with a lower Hashin-Strikmann bound:
# 
# .. math::
#       K_{\mathrm{dry}}=\left[\frac{\phi / \phi_{\mathrm{b}}}{K_{\mathrm{b}}+(4 / 3) \mu_{\mathrm{b}}}+\frac{1-\phi / \phi_{\mathrm{b}}}{K+(4 / 3) \mu_{\mathrm{b}}}\right]^{-1}-\frac{4}{3} \mu_{\mathrm{b}}
# 
# 
# .. math::
#       \mu_{\mathrm{dry}}=\left[\frac{\phi / \phi_{\mathrm{b}}}{\mu_{\mathrm{b}}+z}+\frac{1-\phi / \phi_{\mathrm{b}}}{\mu+z}\right]^{-1}-z, z=\frac{\mu_{\mathrm{b}}}{6}\left(\frac{9 K_{\mathrm{b}}+8 \mu_{\mathrm{b}}}{K_{\mathrm{b}}+2 \mu_{\mathrm{b}}}\right)
# 
#
# Increasing cement model 
# ~~~~~~~~~~~~~~~~~~~~~~~
# The contact cement model represents the initial stage of the “diagenetic trend” in the data. It is found to be applicable to high-porosity sands. During more severe cementation where the diagenetic cement is filling up the pore space, the contact theory breaks down. One should use  the modified Hashin–Shtrikman upper bound (also referred to as the “increasingcement model). The high-porosity end member is determined by contact
# theory. The first 2–3% cement should be modeled with the contact-cement model. Further increase in cement volume and decrease in porosity is described by an HS upper bound interpolation between the high-porosity end member and the mineral point. 
# 
#
# Examples
# ~~~~~~~~
# Let's compute effective bulk and shear moduli of a water saturated rock using different bound models.
#

# %%


# specify model parameters
phic=0.4 # critical porosity
sigma=20 # effective pressure 
Cn=8.6 #coordination number
f=0.5# reduced shear factor
phi = np.linspace(1e-7,phic,100) #define porosity range according to critical porosity
K0, G0 = 36.6, 45 ## grain density, bulk and shear modulus 
Kc, Gc = 36.6, 45 # cement density, bulk and shear modulus
vsh=0 # shale volume
phib=0.3 # critical cement porosity
## softsand, stiffsand and contact cement model
Kdry1, Gdry1 = GM.softsand(K0, G0, phi, phic, Cn, sigma,f) 
Kdry2, Gdry2 = GM.stiffsand(K0, G0, phi, phic, Cn, sigma,f)
Kdry3, Gdry3 = GM.contactcement(K0, G0, Kc, Gc, phi, phic, Cn,  scheme=2)
# plot
plt.figure(figsize=(6,6))
plt.xlabel('Porosity')
plt.ylabel('Bulk modulus [GPa]')
plt.title('')
plt.plot(phi, Kdry1,label='Softsand')
plt.plot(phi, Kdry2,label='Stiffsand')
plt.plot(phi, Kdry3,label='Contact cement')
plt.legend(loc='best')
plt.grid(ls='--')

# %%
# From the figure we can see that the contact cement model will fail for small porosity sand, so we define a critical cement porosity :math:`\phi_b` below which constant cement model/increasing cement model can be applied. 
#

# %%


## constant cement model
Kdry3[phi<phib]=np.nan
Kdry4, Gdry4=GM.constantcement(phib, K0, G0,Kc,Gc, phi, phic, Cn, scheme=2)
Kdry4[phi>phib]=np.nan
## increasing cement model
Kdry5, Gdry5 = GM.MUHS(K0, G0, Kc,Gc,phi, phib,phic, Cn, scheme=2)
Kdry5[phi>phib]=np.nan
#plot
plt.figure(figsize=(6,6))
plt.xlabel('Porosity')
plt.ylabel('Bulk modulus [GPa]')
plt.title('')
plt.plot(phi, Kdry1,label='Softsand')
plt.plot(phi, Kdry2,label='Stiffsand')
plt.plot(phi, Kdry3,label='Contact cement')
plt.plot(phi, Kdry4,label='constant cement')
plt.plot(phi, Kdry5,label='Increasing cement')
plt.legend(loc='best')
plt.grid(ls='--')


# %%
# **References**:
# 
# - Mavko, G., Mukerji, T. and Dvorkin, J., 2020. The rock physics handbook. Cambridge university press.
# 
# - Avseth, P., Mukerji, T. and Mavko, G., 2010. Quantitative seismic interpretation: Applying rock physics tools to reduce interpretation risk. Cambridge university press.
# 
