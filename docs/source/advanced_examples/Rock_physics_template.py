"""
Rock Physics Template (RPT)
===========================
"""
# %%


import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams['font.size']=14
plt.rcParams['font.family']='arial'


# %%


#import rockphypy # import the module 
from rockphypy import GM, Fluid 


# %%
# Rock physics models can be used to calculate elastic properties with various combination of lithology and fluid parameters. Rock Physics Templates (RPTs) were first presented by Ødegaard and Avseth (2003).  Rock physics templates (RPT) is used to display a reference framework of all the possible variations of a particular rock and use such templates to understand actual well log data (or seismic-derived elastic properties). 
# 
# First we write a function to create RPT plot given computed Impedance and Vp/Vs ratio using different rock physics models
#


# %%


def plot_rpt(Kdry,Gdry,K0,D0,Kb,Db,Khc,Dhc,phi,sw):
    """Create RPT plot given computed Impedance and Vp/Vs ratio. 

    Parameters
    ----------
    Kdry : float or array
        effective bulk modulus given by rock physics model
    Gdry : float or array
        effective shear modulus given by rock physics model
    K0 : float
        bulk modulus of grain
    D0 : float
        density of grain
    Kb : float
        bulk moduluf of brine 
    Db : float
        density of brine
    Khc : float
        bulk modulus of HC
    Dhc : float
        density of HC
    phi : float or array 
        porosity
    sw : float or array
        water saturation

    Returns
    -------
    python onject: fig
        rpt plot
    """    
    
    # setup empty arrays to store Ip and Vp/Vs values
    IP=np.empty((phi.size,sw.size))
    PS=np.empty((phi.size,sw.size))

    ## loop over Sw, computes elastic moduli of fluid mixture and saturated rock properties with Gassmann's equation
    
    for i,val in enumerate(sw):
        Kf=(val/Kb+(1-val)/Khc)**-1
        #Kf_mix(val,Kb,Khc)
        Df = val*Db+(1-val)*Dhc
        vp,vs,rho= Fluid.vels(Kdry,Gdry,K0,D0,Kf,Df,phi)
        IP[:,i]=vp*rho
        PS[:,i]=vp/vs
    # plot
    fig=plt.figure(figsize=(10,8))
    plt.plot(IP.T, PS.T, '-ok', mec='k', ms=10, mfc='yellow')
    plt.plot(IP[:,-1], PS[:,-1], '-ok', mec='k', ms=10, mfc='blue')# Brine line 

    plt.xlabel('Acoustic Impedance'), plt.ylabel('Vp/Vs')
    for i,val in enumerate(phi):
        plt.text(IP[i,-1],PS[i,-1]+.03,'{:.02f}'.format(val), alpha=1,backgroundcolor='0.9')
    for i,val in enumerate(sw):
        plt.text(IP[-1,i]-100,PS[-1,i],'Gas={:.02f}'.format(1-sw[i]),ha='right',alpha=1)
    return fig


# %%



# specify model parameters
D0, K0, G0 = 2.65, 36.6, 45
Dc,  Kc, Gc = 2.65,37, 45 # cement
Db, Kb = 1, 2.5
Do, Ko = 0.8, 1.5
Dg, Kg = 0.2, 0.05
### adjustable para
phi_c = 0.4
Cn=8.6  ## calculate coordination number 
phi = np.linspace(0.1,phi_c,10) #define porosity range according to critical porosity
sw=np.linspace(0,1,5) # water saturation
sigma=20
f=0.5

# %%
# - Case 1: create RPT for unconsolidated sand using softsand/friable sand model 
#


# %%


# softsand model gas
Kdry1, Gdry1 = GM.softsand(K0, G0, phi, phi_c, Cn, sigma,f) # soft sand 
fig1=plot_rpt(Kdry1,Gdry1,K0,D0,Kb,Db,Kg,Dg,phi,sw) 
plt.title('Softsand RPT-gas')  
plt.xlim(1000,10000)
plt.ylim(1.4,2.4)


# %%


# softsand model oil
fig1_=plot_rpt(Kdry1,Gdry1,K0,D0,Kb,Db,Ko,Do,phi,sw) 
plt.title('Softsand RPT-oil')  
plt.xlim(1000,10000)
plt.ylim(1.4,2.4)

# %%
# - Case 2: create RPT for stiff sandstone using stiffsand model 
#

# %%


Kdry2, Gdry2 = GM.stiffsand(K0, G0, phi, phi_c, Cn, sigma, f) # stiff sand
fig2=plot_rpt(Kdry2,Gdry2,K0,D0,Kb,Db,Kg,Dg,phi,sw) 
plt.title('Stiffsand RPT-gas')  
plt.xlim(1000,14000)
plt.ylim(1.4,2.3)


# %%


# stiffsand model oil
fig2_=plot_rpt(Kdry2,Gdry2,K0,D0,Kb,Db,Ko,Do,phi,sw) 
plt.title('Stiffsand RPT-oil')  
plt.xlim(1000,14000)
plt.ylim(1.4,2.3)


# %%
# **Reference**: 
#
# - Mavko, G., Mukerji, T. and Dvorkin, J., 2020. The rock physics handbook. Cambridge university press.
#
# - Avseth, P.A. and Odegaard, E., 2003. Well log and seismic data analysis using rock physics templates. First break, 22(10).
#
# - Avseth, P., Mukerji, T. and Mavko, G., 2010. Quantitative seismic interpretation: Applying rock physics tools to reduce interpretation risk. Cambridge university press.
#     

