#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   model.py
@Time    :   2021/06/28 21:20:10
@Author  :   Jiaxin Yu 
@Contact :   yujiaxin666@outlook.com
@License :   (C)Copyright 2020-2021, Jiaxin Yu
'''
import numpy as np
#--------------------------01--------------------------#
def cripor(K0, G0, phi, phic):  
    """Critical porosity model according to Nur’s hypothesis. Written by Jiaxin Yu (July 2021)
    
    Args:
        K0 (GPa): mineral bulk modulus 
        G0 (Gpa): mineral shear modulus 
        phi (frac): porosity
        phic (frac): critical porosity
    Returns:
        K_dry,G_dry: dry elastic moduli of the framework
    """    
    K_dry = K0 * (1-phi/phic)
    G_dry = G0 * (1-phi/phic)

    return K_dry, G_dry

def Gassmann(K_dry,G_dry,K_mat,Kf,phi):
    """Computes velocities and elastic moduli of saturated rock via Gassmann equation given dry-rock moduli. Written by Jiaxin Yu (July 2021)

    Args:
        K_dry (Gpa): dry frame bulk modulus 
        G_dry (Gpa): dry frame shear modulus 
        K_mat (Gpa): matrix bulk modulus
        Kf (Gpa): fluid bulk modulus
        phi (frac): porosity

    Returns:
        K_sat, G_sat: fluid saturated elastic moduli
    """    
    A=(1-K_dry/K_mat)**2
    B=phi/Kf+(1-phi)/K_mat-K_dry/(K_mat**2)
    K_sat=K_dry+A/B
    G_sat = G_dry # At low frequencies, Gassmann’s relations predict no change in the shear modulus between dry and saturated patches
    return K_sat,G_sat
#--------------------------02--------------------------#
def VRH(volumes,M):
    """Computes Voigt-Reuss-Hill Average Moduli Estimate. Written by Jiaxin Yu (July 2021)

    Args:
        volumes (array): volumetric fractions of N phases
        M (array): elastic modulus of the N phase.
    Returns:
        M_v: Voigt average
        M_r: Reuss average
        M_0: Hill average 
    """    
    M_v=np.dot(volumes,M)
    
    M_r=np.dot(volumes,1/M)**-1
    
    M_h= 0.5*(M_r+M_v)
    return  M_v,M_r,M_h

def HS(f, K1, K2,G1, G2, bound='upper'): 
    """Compute two-phase hashin-strikmann bound. Written by Jiaxin Yu (July 2021)

    Args:
        f (float): 0-1, volume fraction of stiff material  
        K1 (GPa): bulk modulus of stiff phase
        K2 (GPa): bulk modulus of soft phase
        G1 (GPa): shear modulus of stiff phase
        G2 (GPa): shear modulus of soft phase
        bound (str, optional): upper bound or lower bound. Defaults to 'upper'.
    """    
    if bound == 'upper':  
        K=K1+ (1-f)/( (K2-K1)**-1 + f*(K1+4*G1/3)**-1 )

        Temp= (K1+2*G1)/(5*G1 *(K1+4*G1/3))
        G=G1+(1-f)/( (G2-G1)**-1 + 2*f*Temp)
    else:  
        K=K2+ f/( (K1-K2)**-1 + (1-f)*(K2+4*G2/3)**-1 )

        Temp= (K2+2*G2)/(5*G2 *(K2+4*G2/3))
        G=G2+f/( (G1-G2)**-1 + 2*(1-f)*Temp)
    return K, G
#--------------------------03--------------------------#
def Dilute_incl(Ks,Gs,phi):
    """Compute effective elastic moduli via "Swiss cheese" model with spherical pores. Written by Jiaxin Yu (July 2021)

    Args:
        Ks (GPa): Bulk modulus of matrix 
        Gs (GPa): shear modulus of matrix
        phi (frac): porosity

    Returns:
        Kdry,Gdry (GPa): effective elastic moduli
    """    
    Kdry=(1/Ks *(1+(1+3*Ks/(4*Gs))*phi))**-1
    Gdry=(1/Gs * (1+(15*Ks+20*Gs)*phi/(9*Ks+8*Gs)))**-1
    return Kdry, Gdry

def SC(phi,Ks,Gs,iter_n):
    """
    Self-Consistent(SC) EM model with spherical pores. Written by Jiaxin Yu (July 2021)

    Args:
        phi (float): volumetric fraction
        Ks (GPa): bulk modulus of matrix phase
        Gs (GPa): shear modulus of matrix phase
        
        iter_n (int): iterations, necessary iterations increases as f increases.
    Note: 
        phi.shape== Km.shape
    """  
    K_eff=Ks
    G_eff=Gs
    for i in range(iter_n):
        K_eff = (1/Ks + (1/K_eff+3/(4*G_eff))*phi) **-1
        G_eff= (1/Gs + (15*K_eff+20*G_eff)/(9*K_eff+8*G_eff) *phi/G_eff )**-1
    return K_eff,G_eff

#--------------------------04--------------------------#
def hertzmindlin(K0, G0, phic, Cn, sigma, f):
    """Compute effective elastic moduli of granular packing via Hertz-Mindlin approach. Written by Jiaxin Yu (July 2021)
    Args:
        K0 (GPa): bulk modulus of grain material
        G0 (GPa): shear modulus of grain material
        phic (GPa): critical porosity
        Cn (unitless): coordination number
        sigma (MPa): effective stress
        f (unitless): reduced shear factor, 0=dry pack with inifinitely rough spheres; 1=dry pack with infinitely smooth spheres

    Returns:
        K_dry, G_dry: effective elastic moduli of dry pack
    
    """    
    sigma =sigma/1000 # converts pressure unit to GPa
    nu=(3*K0-2*G0)/(6*K0+2*G0) # poisson's ratio of mineral mixture
    K_dry = (sigma*(Cn**2*(1-phic)**2*G0**2) / (18*np.pi**2*(1-nu)**2))**(1/3)
    G_dry = ((2+3*f-nu*(1+3*f))/(5*(2-nu))) * ((sigma*(3*Cn**2*(1-phic)**2*G0**2)/(2*np.pi**2*(1-nu)**2)))**(1/3)
    return K_dry, G_dry

def Walton(K0, G0, phic, Cn, sigma, f):
    """ Compute dry rock elastic moduli within the mechanical compaction
    domain based on the Contact Theory (CT) of Walton (1987). Written by Jiaxin Yu (July 2021)

    Args:
        K0 (GPa): Bulk modulus of grain material
        G0 (GPa): shear modulus of grain material
        phic (frac): depositional porosity, Defaults to 0.4.
        Cn (unitless): coordination number Defaults to 8.6.
        sigma (MPa): effective stress. Defaults to 10.
        f (unitless): reduced shear factor. 0=dry pack with inifinitely rough spheres; 1=dry pack with infinitely smooth spheres
    """   
    sigma= sigma/1e3 # convert Mpa to Gpa 
    lamda = K0- 2*G0/3  # Lamé’s coefficient of the grain material.  
    B = 1/(4*np.pi) *(1/G0+1/(G0+lamda))
    #A = 1/(4*np.pi) *(1/G0-1/(G0+lamda))

    K_w= 1/6 * (3*(1-phic)**2*Cn**2*sigma/(np.pi**4*B**2))**(1/3)
    G_dry1= 3/5 * K_w
    nu =(3*K0-2*G0)/(6*K0+2*G0)  # possion's ratio
    G_dry2= 3/5 * K_w *(5-4*nu)/(2-nu) # rought limit is the same as Hertz-Mindlin
    G_w= f* G_dry2+(1-f)*G_dry1
    return K_w, G_w

#--------------------------06--------------------------#
def M_from_log(den, vp,vs):
    '''
    Compute K_sat  and G from well log data. Written by Jiaxin Yu (July 2021)
    Parameters:
    den: bulk density, DEN
    vp,vs
    
    '''
    K=den*1000*(vp**2-4*vs**2/3)
    G=den*1000*vs**2
    K=K/10**9
    G=G/10**9
    return K, G

def Thomsen(C11, C33, C13, C44, C66, den, theta):
    """Computes  three phase velocities for weak anisotropy using Thomsen’s parameters. Written by Jiaxin Yu (July 2021)

    Args:
        C11 (GPa): stiffness
        C33 (GPa): stiffness
        C13 (GPa): stiffness
        C44 (GPa): stiffness
        C66 (GPa): stiffness
        den (g/cm3): density of the medium
        theta (unitless): angle of incidence

    Returns:
        VP, VSV, VSH: wave velocities propagating along given direction
    """    
    alpha= np.sqrt(C33/den) *1e3
    beta= np.sqrt(C44/den) *1e3
    epsilon= 0.5* (C11-C33)/C33
    gamma= 0.5 * (C66-C44)/C44
    theta= np.deg2rad(theta)
    delta= ((C13+C44)**2-(C33-C44)**2 )/ (2*C33*(C33-C44))
    VP= alpha*(1+delta*np.sin(theta)**2*np.cos(theta)**2+epsilon*np.sin(theta)**4)
    VSV= beta*(1+alpha**2/beta**2 * (epsilon-delta)*np.sin(theta)**2*np.cos(theta)**2)
    VSH= beta*(1+gamma*np.sin(theta)**2)
    return VP, VSV, VSH

def Backus(V,lamda, G ):
    """Computes stiffnesses of a layered medium using backus average model. Written by Jiaxin Yu (July 2021)

    Args:
        V (num or array-like, frac): volumetric fractions of N isotropic layering materials
        lamda (num or array-like): Lamé coefficients of N isotropic layering materials
        G (num or array-like, GPa): shear moduli of N isotropic layering materials
    Returns:
        C11,C33,C13,C44,C66 (num or array-like, GPa): Elastic moduli of the anisotropic layered media
    """    
    C33=np.dot(V, 1/(lamda+2*G)) **-1
    C44=np.dot(V, 1/G)**-1
    C66=np.dot(V, G)
    C13=np.dot(V, 1/(lamda+2*G)) **-1 * np.dot(V, lamda/(lamda+2*G))
    C11=np.dot(V, 4*G*(lamda+G)/(lamda+2*G))+np.dot(V, 1/(lamda+2*G))**-1 * np.dot(V, lamda/(lamda+2*G))**2
    
    return C11,C33,C13,C44,C66
#--------------------------07--------------------------#
def exp_decay(c, phic, Z):
    """compute porosity reduction caused by mechanical compaction using Arthy equation. Written by Jiaxin Yu (July 2021)

    Args:
        c (1/m): exponentional decay parameter
        phic (frac): depositional porosity
        Z (m): burial depth

    Returns:
        phi (frac): porosity
    """    
    phi= phic*np.exp(-c*Z)
    return phi

def IGV(IGV_f,phic,beta,sigma ):
    """model porosity evolution of a clean, uncemented sand caused by mechanical compaction using the Lander and Walderhaug (1999) equations. Written by Jiaxin Yu (July 2021)

    Args:
        IGV_f (frac): stable packing configuration, which is the smallest porosity possible by mechanical compaction alone.
        phic (frac): depositional porosity 
        beta (1/MPa): exponential rate of IGV decline with effective stress, the higher, the more severe mechanical compaction the rock experienced. 0.06 is taken in Torset （2020) and 0.12  is taken in Bredesen（2017).
        sigma (MPa):  maximum effective stress 
    Returns:
        IGV (frac): Intergranular volume: part of the rock between grains
    """    
    IGV=IGV_f+(phic-IGV_f)*np.exp(-beta*sigma)
    return IGV
def Qz_cement(n,tcem,phic,cc, d=70,D=0.03,f_q=0.85,a=1.981*10**-22,b=0.022,M=60.9,coat=0,den_q=2.65):
    """Quantify volume of quartz precipite in the pore space in a period of time using Walderhaug (1994a, 1996) kinetic equations. Written by Jiaxin Yu (July 2021)

    Args:
        n (int): additive time step, default 1 for 1 m.y timestep
        tcem (m.y.): Cementation time
        phic (frac):  porosity at the start of cementation
        cc: (°C/m.y.):  heating rate. 
        d (°C): Temperature at the onset of cementation. Defaults to 70.
        D (mm): grain diameter. Defaults to 0.03.
        f_q (frac): volume fraction of intergranular quartz in a unit volume V . Defaults to 0.85.
        a (moles/cm2 s): kinetic constants. Defaults to 1.981*10**-22.
        b (°C−1): kinetic constants. Defaults to 0.022.
        M ( g/mole): molar mass of quartz. Defaults to 60.9.
        coat (int, optional):. Grain coatings prohibiting further quartz overgrowth are compensated 
        by the factor. Defaults to 0.
        den_q (g/cm3):  density of quartz. Defaults to 2.65.
    Returns:
        V (frac): volume fraction of quartz cementation
        Por (frac): remaining porosity after cementation
        A (cm2/cm3): surface area remained
    """    
    convert=365*24*60*60*1000000
    a=a*convert # moles/ (cm2 my.)
    A0=6*(1-coat)*f_q/D # V is 1 cm3 initial quartz grains surface area
    tt=np.linspace(0,tcem,(int(tcem)+1)*n) # time step, usually set to 1 my 
    V=np.zeros(len(tt))
    A=np.zeros(len(tt))
    Por=np.zeros(len(tt))
    A[0]=A0
    Por[0]=phic
    for i in range(1,len(tt)):
    
        Term1=-M*a*A0/(den_q*phic*b*cc*np.log(10))
        Term2=10**(b*(cc*tt[i]+d))-10**(b*(cc*tt[i-1]+d)) 
        V[i]=phic-(phic-V[i-1])*np.exp(Term1*Term2)

        A[i]=A[0]*(1-V[i]/phic)
        Por[i]=phic-V[i]
    return V, Por, A,tt
#--------------------------08--------------------------#
def softsand(K0, G0, phi, phic, Cn, sigma, f):
    """Soft-sand (unconsolidated sand) model: model the porosity-sorting effects using the lower Hashin–Shtrikman–Walpole bound. (Also referred to as the ‘friable-sand model’ in Avseth et al. 2010). Written by Jiaxin Yu (July 2021)

    Args:
        K0 (GPa): Bulk modulus of grain material
        G0 (GPa): shear modulus of grain material
        phi (float or array-like): Porosity
        phic (frac): depositional porosity, Defaults to 0.4.
        Cn (unitless): coordination number Defaults to 8.6.
        sigma (MPa): effective stress. 
        f (unitless): reduced shear factor. 0=dry pack with inifinitely rough spheres; 1=dry pack with infinitely smooth spheres
    Returns:
        K_dry, G_dry (GPa): Effective elastic moduli of dry pack
    """   
    K_HM, G_HM = hertzmindlin(K0, G0, phic, Cn, sigma, f)
    K_dry =-4/3*G_HM + (((phi/phic)/(K_HM+4/3*G_HM)) + ((1-phi/phic)/(K0+4/3*G_HM)))**-1
    aux = G_HM/6*((9*K_HM+8*G_HM) / (K_HM+2*G_HM)) # auxiliary variable 
    G_dry = -aux + ((phi/phic)/(G_HM+aux) + ((1-phi/phic)/(G0+aux)))**-1
    return K_dry, G_dry
def stiffsand(K0, G0, phi, phic, Cn, sigma, f):
    """Stiff-sand model: model the porosity-sorting effects using the lower Hashin–Shtrikman–Walpole bound. Written by Jiaxin Yu (July 2021)

    Args:
        K0 (GPa): Bulk modulus of grain material
        G0 (GPa): shear modulus of grain material
        phi (float or array-like): Porosity
        phic (frac): depositional porosity, Defaults to 0.4.
        Cn (unitless): coordination number Defaults to 8.6.
        sigma (MPa): effective stress. 
        f (unitless): reduced shear factor. 0=dry pack with inifinitely rough spheres; 1=dry pack with infinitely smooth spheres

    Returns:
        K_dry, G_dry (GPa): Effective elastic moduli of dry pack
    """    
    K_HM, G_HM = hertzmindlin(K0, G0, phic, Cn, sigma, f)
    K_dry = -4/3*G0 + (((phi/phic)/(K_HM+4/3*G0)) + ((1-phi/phic)/(K0+4/3*G0)))**-1
    aux = G0/6*((9*K0+8*G0) / (K0+2*G0))
    G_dry = -aux + ((phi/phic)/(G_HM+aux) + ((1-phi/phic)/(G0+aux)))**-1
    return K_dry, G_dry
def contactcement(K0, G0, Kc, Gc, phi, phic, Cn,  scheme):
    """Compute dry elastic moduli of cemented sandstone via Contact cement (cemented sand) model, Dvorkin-Nur (1996). Written by Jiaxin Yu (July 2021)

    Args:
        K0 (GPa): Bulk modulus of grain material
        G0 (GPa): shear modulus of grain material
        Kc (GPa): Bulk modulus of cement
        Gc (GPa): shear modulus of cement
        phi (float or array-like): Porosity
        phic (frac): depositional porosity, Defaults to 0.4.
        Cn (unitless): coordination number Defaults to 8.6.
        scheme: Scheme of cement deposition
                1=cement deposited at grain contacts
                2=cement deposited at grain surfaces
    Returns:
        K_dry, G_dry (GPa): Effective elastic moduli of dry rock
    """    
    nu_0=(3*K0-2*G0)/(6*K0+2*G0) # Poisson's ratio of grain material
    nu_c = (3*Kc-2*Gc)/(6*Kc+2*Gc) # Poisson's ratio of cement
    if scheme == 1: # scheme 1: cement deposited at grain contacts
        alpha = 2*((phic-phi)/(3*Cn*(1-phic))) ** (1/4)
    else: # scheme 2: cement evenly deposited on grain surface
        alpha = ((2*(phic-phi))/(3*(1-phic)))**(1/2)
    LambdaN = (2*Gc*(1-nu_0)*(1-nu_c)) / (np.pi*G0*(1-2*nu_c))
    N1 = -0.024153*LambdaN**-1.3646
    N2 = 0.20405*LambdaN**-0.89008
    N3 = 0.00024649*LambdaN**-1.9864
    Sn = N1*alpha**2 + N2*alpha + N3
    LambdaT = Gc/(np.pi*G0)
    T1 = -10**-2*(2.26*nu_0**2+2.07*nu_0+2.3)*LambdaT**(0.079*nu_0**2+0.1754*nu_0-1.342)
    T2 = (0.0573*nu_0**2+0.0937*nu_0+0.202)*LambdaT**(0.0274*nu_0**2+0.0529*nu_0-0.8765)
    T3 = 10**-4*(9.654*nu_0**2+4.945*nu_0+3.1)*LambdaT**(0.01867*nu_0**2+0.4011*nu_0-1.8186)
    St = T1*alpha**2 + T2*alpha + T3
    K_dry = 1/6*Cn*(1-phic)*(Kc+(4/3)*Gc)*Sn
    G_dry = 3/5*K_dry+3/20*Cn*(1-phic)*Gc*St
    return K_dry, G_dry
def constantcement(phi_b, K0, G0, Kc, Gc,phi, phic, Cn,scheme):
    """Constant cement (constant depth) model according to Avseth (2000)
       https://doi.org/10.1190/1.3483770. Written by Jiaxin Yu (July 2021)

    Args:
        phi_b (frac): adjusted high porosity end memeber
        K0 (GPa): Bulk modulus of grain material
        G0 (GPa): shear modulus of grain material
        Kc (GPa): Bulk modulus of cement
        Gc (GPa): shear modulus of cement
        phi (float or array-like): Porosity
        phic (frac): depositional porosity, Defaults to 0.4.
        Cn (unitless): coordination number Defaults to 8.6.
        scheme: Scheme of cement deposition
                1=cement deposited at grain contacts
                2=cement deposited at grain surfaces
    Returns:
        K_dry, G_dry (GPa): Effective elastic moduli of dry rock
    """    
    K_b,G_b=contactcement(K0, G0, Kc, Gc, phi_b, phic, Cn,  scheme)
    T=phi/phi_b
    Z=(G_b/6)*(9*K_b+8*G_b)/(K_b+2*G_b)
    Kdry=(T/(K_b+4*G_b/3)+(1-T)/(K0+4*G_b/3))**(-1)-4*G_b/3
    Gdry=(T/(G_b+Z)+(1-T)/(G0+Z))**(-1)-Z

    return Kdry, Gdry  
    
def MUHS(K0, G0, Kc,Gc,phi, phi_b,phic, Cn,scheme):
   """Increasing cement model: Modified Hashin-Strikmann upper bound blend with contact cement model. Written by Jiaxin Yu (July 2021)
   
   Args:
      K0 (GPa): Bulk modulus of grain material
      G0 (GPa): shear modulus of grain material
      Kc (GPa): Bulk modulus of cement
      Gc (GPa): shear modulus of cement
      phi (float or array-like): Porosity
      phi_b (frac): adjusted high porosity end memeber
      phic (frac): depositional porosity, Defaults to 0.4.
      Cn (unitless): coordination number Defaults to 8.6.
      scheme: Scheme of cement deposition
                1=cement deposited at grain contacts
                2=cement deposited at grain surfaces
   Returns:
      K_dry, G_dry (GPa): Effective elastic moduli of dry rock
   """   
   # adjusted high porosity end memeber 
   K_b, G_b=contactcement(K0, G0, Kc, Gc, phi_b, phic, Cn, scheme)
   # interpolation between the high-porosity end member Kb, Gb and the mineral point. 
   K_dry = -4/3*G0 + (((phi/phi_b)/(K_b+4/3*G0)) + ((1-phi/phi_b)/(K0+4/3*G0)))**-1
   tmp = G0/6*((9*K0+8*G0) / (K0+2*G0))
   G_dry = -tmp + ((phi/phi_b)/(G_b+tmp) + ((1-phi/phi_b)/(G0+tmp)))**-1
   return K_dry, G_dry

#--------------------------09--------------------------#
def silty_shale(C, Kq,Gq, Ksh, Gsh):
    """Dvorkin–Gutierrez silty shale model: model the elastic moduli of decreasing clay content for shale. Written by Jiaxin Yu (July 2021)

    Args:
        C (frac): volume fraction of clay
        Kq (GPa): bulk modulus of silt grains
        Gq (GPa): shear modulus of silt grains
        Ksh (GPa): saturated bulk modulus of pure shale
        Gsh (GPa): saturated shear modulus of pure shale
        ** Ksh and Gsh could be derived from well-log measurements of VP, VS and density in a pure shale zone.

    Returns:
        K_sat, G_sat: elastic moduli of the saturated silty shale.
    """    

    K_sat = ( C/(Ksh+4*Gsh/3)+ (1-C)/(Kq+4*Gq/3) )**-1 -4*Gsh/3
    Zsh = Gsh/6 *(9*Ksh+8*Gsh)/(Ksh+2*Gsh)
    G_sat = (C/(Gsh+Zsh) + (1-C)/(Gq+Zsh))**-1 -Zsh
    return K_sat,G_sat

def shaly_sand(phis, C, Kss,Gss, Kcc, Gcc):
    """Modeling elastic moduli for sand with increasing clay content using LHS bound rather than using Gassmann relation. Written by Jiaxin Yu (July 2021)

    Args:
        phis (porosity): critical porosity of sand composite
        C (frac): clay content 
        Kss (GPa): saturated bulk moduli for clean sandstone using e.g. HM
        Gss (GPa): saturated shear moduli for clean sandstone using e.g. HM
        Kcc (GPa): saturated bulk moduli calculated from the sandy shale model at critical clay content using silty shale model  
        Gcc (GPa): saturated shear moduli calculated from the sandy shale model at critical clay content using silty shale model

    Returns:
        K_sat,G_sat: saturated rock moduli of the shaly sand 
    """    
    K_sat = ( (1-C/phis)/(Kss+4*Gss/3)+ (C/phis)/(Kcc+4*Gss/3) )**-1 -4*Gss/3
    Zss = Gss/6 *(9*Kss+8*Gss)/(Kss+2*Gss)
    G_sat = ((1-C/phis)/(Gss+Zss) + (C/phis)/(Gcc+Zss))**-1 -Zss
    return K_sat,G_sat

#--------------------------10--------------------------#
def vels(K_dry,G_dry,K0,D0,Kf,Df,phi):
    """ Computes Vp,Vs and densities of saturated rock using Gassmann relations. Written by Jiaxin Yu (July 2021)

    Args:
        K_dry (GPa): dry frame bulk modulus
        G_dry (GPa): dry frame shear modulus
        K0 (GPa): mineral matrix bulk modulus
        D0 (GPa): mineral matrix density
        Kf (GPa): fluid bulk modulus
        Df (g/cm3): fluid density
        phi (frac): porosity

    Returns:
        [type]: [description]
    """    
    #D0=Dsh*vsh+Dqz*(1-vsh)####
    rho  = D0*(1-phi)+Df*phi
    K,_= Gassmann(K_dry,G_dry,K0,Kf,phi)
    Vp   = np.sqrt((K+4./3*G_dry)/rho)*1e3
    Vs   = np.sqrt(G_dry/rho)*1e3
    return Vp, Vs, rho
def plot_rpt(Kdry,Gdry,K0,D0,Kb,Db,Khc,Dhc,phi,sw):
    """Create RPT plot given computed Impedance and Vp/Vs ratio. Written by Jiaxin Yu (July 2021)

    Args:
        IP (2D array): Impedance
        PS (2D array): Vp/Vs ratio
        phi (): [description]
        sw ([type]): [description]

    Returns:
        [type]: [description]
    """    
    # setup empty arrays to store Ip and Vp/Vs values
    IP=np.empty((phi.size,sw.size))
    PS=np.empty((phi.size,sw.size))

    ## loop over Sw, computes elastic moduli of fluid mixture and saturated rock properties with Gassmann's equation
    #(Khc, Dhc) = (Kg, Dg) if fluid == 'gas' else (Ko,Do)
    for i,val in enumerate(sw):
        Kf=(val/Kb+(1-val)/Khc)**-1
        #Kf_mix(val,Kb,Khc)
        Df = val*Db+(1-val)*Dhc
        vp,vs,rho= vels(Kdry,Gdry,K0,D0,Kf,Df,phi)
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

#--------------------------11--------------------------#
from scipy.integrate import odeint

def PQ(Km,Gm, Ki,Gi, alpha):    
    """ compute geometric strain concentration factors P and Q for spheroids of arbitrary aspect ratio according to Berymann (1980). Written by Jiaxin Yu (July 2021)

    Args:
        Km (GPa): Bulk modulus of matrix phase
        Gm (GPa): Shear modulus of matrix phase
        Ki (GPa): Bulk modulus of matrix phase
        Gi (GPa): Shear modulus of matrix phase
        alpha (float): aspect ratio

    Returns:
        P,Q (unitless): geometric strain concentration factors
    """
    theta= alpha/(1.0 - alpha**2)**(3.0/2.0) * (np.arccos(alpha) - alpha*np.sqrt(1.0 - alpha**2))
    f= alpha**2*(3.0*theta - 2.0)/(1.0 - alpha**2)
    A = Gi/Gm - 1.0
    B = (Ki/Km - Gi/Gm)/3.0
    R = Gm/(Km + (4.0/3.0)*Gm) # 
    F1 = 1.0 + A*(1.5*(f + theta) - R*(1.5*f + 2.5*theta - 4.0/3.0))
    F2 = 1.0 + A*(1.0 + 1.5*(f + theta) - R*(1.5*f + 2.5*theta)) + B*(3.0 - 4.0*R) + A*(A + 3.0*B)*(1.5 - 2.0*R)*(f + theta - R*(f - theta + 
    2.0*theta**2))
    F3 = 1.0 + A*(1.0 - f - 1.5*theta + R*(f + theta))
    F4 = 1.0 + (A/4.0)*(f + 3.0*theta - R*(f - theta))
    F5 = A*(-f + R*(f + theta - 4.0/3.0)) + B*theta*(3.0 - 4.0*R)
    F6 = 1.0 + A*(1.0 + f - R*(f + theta)) + B*(1.0 - theta)*(3.0 - 4.0*R)
    F7 = 2.0 + (A/4.0)*(3.0*f + 9.0*theta - R*(3.0*f + 5.0*theta)) + B*theta*(3.0 - 4.0*R)
    F8 = A*(1.0 - 2.0*R + (f/2.0)*(R - 1.0) + (theta/2.0)*(5.0*R - 3.0)) + B*(1.0 - theta)*(3.0 - 4.0*R)
    F9 = A*((R - 1.0)*f - R*theta) + B*theta*(3.0 - 4.0*R)
    Tiijj = 3*F1/F2
    Tijij = Tiijj/3 + 2/F3 + 1/F4 + (F4*F5 + F6*F7 - F8*F9)/(F2*F4)
    P = Tiijj/3
    Q = (Tijij - P)/5
    return P, Q
def DEM(y,t, params):
    '''
    ODE solver tutorial: https://physics.nyu.edu/pine/pymanual/html/chap9/chap9_scipy.html. Written by Jiaxin Yu (July 2021)
    '''
    K_eff,G_eff=y  # unpack current values of y
    Gi,Ki,alpha = params # unpack parameters 
    P, Q= PQ(G_eff,K_eff,Gi,Ki, alpha)
    derivs = [1/(1-t) * (Ki-K_eff) * P,  1/(1-t) * -G_eff * Q]
    return derivs

#--------------------------12--------------------------#
def PP_ref(theta, vp1,vp2,vs1,vs2,den1,den2):
    """Aki-Richard approximation to PP reflectivity. Written by Jiaxin Yu (July 2021)

    Args:
        theta (degree): incident angle
        vp1 (m/s): P wave velocity of layer 1
        vp2 (m/s): P wave velocity of layer 2
        vs1 (m/s): S wave velocity of layer 1
        vs2 (m/s): S wave velocity of layer 2
        den1 (kg/m3): density of layer 1
        den2 (kg/m3): density of layer 2

    Returns:
        R_pp: P wave reflectivity
    """   
    theta=np.deg2rad(theta) # convert angle in degree to angle in radian
    delta_den=den2-den1
    delta_vp=vp2-vp1
    delta_vs=vs2-vs1
    rho_mean=0.5*(den1+den2)
    vp_mean=0.5*(vp1+vp2)
    vs_mean=0.5*(vs1+vs2)
    
    Rpp0=0.5*(delta_den/rho_mean+delta_vp/vp_mean)
    M= -2*(vs_mean/vp_mean)**2*(2*delta_vs/vs_mean+delta_den/rho_mean)
    N= 0.5* delta_vp/vp_mean
    R_pp= Rpp0+ M *np.sin(theta)**2 + N * np.tan(theta)**2
    return R_pp

def R_G(vp1,vp2,vs1,vs2,den1,den2):
    """Compute gradient and intercept based on two term approximation. Written by Jiaxin Yu (July 2021)
    Args:
        vp1 (m/s): P wave velocity of layer 1
        vp2 (m/s): P wave velocity of layer 2
        vs1 (m/s): S wave velocity of layer 1
        vs2 (m/s): S wave velocity of layer 2
        den1 (kg/m3): density of layer 1
        den2 (kg/m3): density of layer 2

    Returns:
        R: intercept
        G: gradient
    """  
    delta_den=den2-den1
    delta_vp=vp2-vp1
    delta_vs=vs2-vs1
    rho_mean=0.5*(den1+den2)
    vp_mean=0.5*(vp1+vp2)
    vs_mean=0.5*(vs1+vs2)
    
    R=0.5*(delta_den/rho_mean+delta_vp/vp_mean)
    G= 0.5*delta_vp/vp_mean-2*(vs_mean/vp_mean)**2*(2*delta_vs/vs_mean+delta_den/rho_mean)
    
    return R,G
#--------------------------13--------------------------#
