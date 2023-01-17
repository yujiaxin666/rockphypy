#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   utils.py
@Time    :   2023/01/16 12:06:24
@Author  :   Jiaxin Yu 
@Contact :   yujiaxin666@outlook.com
@License :   (C)Copyright 2020-2021, Jiaxin Yu
'''
# moduli transfor and velocity calculation

import numpy as np


def V(K, G, rho):
    """Compute velocity given density and elastic moduli. 
    Args:
        K (GPa): bulk modulus
        G (GPa): shear moulus
        rho (g/m3): density of the frame 

    Returns:
        Vp, Vs (m/s): velocity
    Written by Jiaxin Yu (July 2021)

    """    
    Vp   = np.sqrt((K+4./3*G)/rho)*1e3
    Vs   = np.sqrt(G/rho)*1e3
    return Vp, Vs

def poi(K, G):
    """Compute poisson's ratio from K an G

    Args:
        K (GPa): bulk modulus
        G (GPa): shear modulus

    Returns:
        nu: Poisson's ratio
    Written by Jiaxin Yu (July 2021)

    """    
    nu=(3*K-2*G)/(6*K+2*G)
    return nu

def lame(K, G):
    """Compute lame constant lamdba from K an G

    Args:
        K (GPa): bulk modulus
        G (GPa): shear modulus

    Returns:
        nu: Poisson's ratio
    Written by Jiaxin Yu (July 2021)

    """    
    lamda= K-2*G/3
    return lamda

def M_from_V(den, vp,vs):
    """Compute K and G from velocities and density

    Args:
        den (g/cm3: bulk density
        vp (m/s): p wave velocity
        vs (m/s): s wave velocity

    Returns:
        K, G (GPa):bulk and shear moduli
    Written by Jiaxin Yu (July 2021)
    """    

    K=den*1000*(vp**2-4*vs**2/3)
    G=den*1000*vs**2
    K=K/10**9
    G=G/10**9
    return K, G

def write_HTI_matrix(C11,C33,C13,C44,C55):
    """formulate HTI stiffness matrix 

    Args:
        C11 (GPa): stiffness
        C13 (GPa): stiffness
        C23 (GPa): stiffness
        C33 (GPa): stiffness
        C44 (GPa): stiffness
        C55 (GPa): stiffness

    Returns:
        C: 6x6 stiffness matrix
    """   
    C23= C33-2*C44
    C=np.array([[C11,C13,C13,0,0,0],
            [C13,C33,C23,0,0,0],
            [C13,C23,C33,0,0,0],
            [0,0,0,C44,0,0],
            [0,0,0,0,C55,0],
            [0,0,0,0,0,C55]])
    return C

def write_VTI_matrix(C11,C33,C13,C44,C66):
    """formulate VTI stiffness matrix 

    Args:
        C11 (GPa): stiffness
        C33 (GPa): stiffness
        C13 (GPa): stiffness
        C44 (GPa): stiffness
        C65 (GPa): stiffness

    Returns:
        C: 6x6 stiffness matrix
    """   
    C12= C11-2*C66
    C=np.array([[C11,C12,C13,0,0,0],
            [C12,C11,C13,0,0,0],
            [C13,C13,C33,0,0,0],
            [0,0,0,C44,0,0],
            [0,0,0,0,C44,0],
            [0,0,0,0,0,C66]])
    return C



