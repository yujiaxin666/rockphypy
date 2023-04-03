"""
Rock physics data screening
===========================
"""

# %%


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['font.size']=14
plt.rcParams['font.family']='arial'
plt.rcParams['axes.labelpad'] = 10.0


# %%


# import the module 
from rockphypy import QI


# %%


import matplotlib.colors
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","orange","yellow"])

# %%
# Given tons of petrophysical properties meaured or derived e.g. from well logs, data screening is inevitably required for later intepretation and scenerial modelling. The purpose of data screening is to identify and address any errors, inconsistencies, or missing data in the dataset before any analysis is conducted. 
#
# Avseth et al. (2020) proposed a rock physics model-based data screening approach. Several contact based elastic models are used to generate elastic bounds for high-porosity sands and sandstones (typical reservoir sandstone). These elastic bounds can tell if data comply with physics. Rock physics diagnostic models used are ``friable sand model``,``contact cement model`` and ``increasing cement model``, ``rockphypy`` provides all the implementations. In addition, ``rockphypy`` has written the code used to plot the elastic bounds as a function ``QI.screening``, ``QI`` stands for Quantitative Intepretation. There are many practically useful functionalities in this submodule.  
#

# %%


# parameters 
Dqz, Kqz, Gqz = 2.65, 36.6, 45 ## grain density, bulk and shear modulus 
Dsh, Ksh, Gsh = 2.7, 21, 7 # shale/clay density, bulk and shear modulus
Dc,Kc, Gc =2.65, 36.6, 45 # cement density, bulk and shear modulus
Db, Kb = 1, 2.2 # brine density, bulk modulus
Do, Ko = 0.8, 1.5 # oil density, bulk modulus
Dg, Kg = 0.2, 0.06 # gas density, bulk modulus
phi_c=0.4 # critical porosity
sigma=20 # effective pressure 
scheme=2
Cn=8.6
# could be array
vsh=0 # shale volume
#phib_p=[0.3,0.36,0.38,0.39] # define cement porosity for Vp
phib_p=0.3
f= 0.5 # slip factor 


# %%


phi,vp1,vp2,vp3,vs1,vs2,vs3 = QI.screening(Dqz,Kqz,Gqz,Dsh,Ksh,Gsh,Dc,Kc,Gc,Db,Kb,phib_p,phi_c,sigma,vsh,scheme,f, Cn)

fig,ax=plt.subplots()
fig.set_size_inches(7, 6)
ax.plot(phi,vp3,'-k', lw=4, alpha=0.7)
ax.plot(phi,vp1,'--k', lw=2, alpha=0.7)
ax.plot(phi,vp2,'-k',lw=4, alpha=0.7)
ax.set_ylabel('Vp (m/s)')
ax.set_xlabel('porosity')
ax.grid(ls='--',alpha=0.7)

# %%
# Applied to field data 
# ^^^^^^^^^^^^^^^^^^^^^
# Let's import a example synthetic well log data and apply the rock physics screening to the well log data 
# 

# %%


data = pd.read_csv('../../data/well/example_well.csv')


# %%

# sphinx_gallery_thumbnail_number = 2
fig,ax=plt.subplots()
fig.set_size_inches(7, 6)
ax.plot(phi,vp3,'-k', lw=4, alpha=0.7)
ax.plot(phi,vp1,'--k', lw=2, alpha=0.7)
ax.plot(phi,vp2,'-k',lw=4, alpha=0.7)
ax.set_ylabel('Vp (m/s)')
ax.set_xlabel('Porosity')
ax.grid(ls='--',alpha=0.7)


plt.scatter(data.PHIT_D,data.VP*1000,c=1-data.VSH_GR,vmin=0, vmax=1,edgecolors='grey',s=100,alpha=1,cmap=cmap1)

cbar=plt.colorbar()
cbar.set_label(r'$V_{\rm Sand}$')

# %%
# **Reference** 
# - Avseth, P., Lehocki, I., Kjøsnes, Ø., & Sandstad, O. (2021). Data‐driven rock physics analysis of North Sea tertiary reservoir sands. Geophysical Prospecting, 69(3), 608-621.
# 