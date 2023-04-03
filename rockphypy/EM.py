#!/usr/bin/env python
# -*-coding:utf-8 -*-


import numpy as np
from rockphypy.utils import utils
#import Anisotropy
from scipy.optimize import fsolve
from scipy.integrate import odeint

class EM:
    """classical bounds and inclusion models. 
    """
    @staticmethod
    def VRH(volumes,M):
        """Computes Voigt, Reuss, and Hill Average Moduli Estimate.

        Parameters
        ----------
        volumes : list or array-like
            volumetric fractions of N phases
        M : list or array-like
            elastic modulus of the N phase.

        Returns
        -------
        float
            M_v: Voigt average
            M_r: Reuss average
            M_0: Hill average
        """        
  
        volumes= np.asanyarray(volumes)
        M=np.asanyarray(M)
        M_v=np.dot(volumes,M)

        M_r=np.dot(volumes,1/M)**-1

        M_h= 0.5*(M_r+M_v)
        return  M_v,M_r,M_h

    @staticmethod
    def cripor(K0, G0, phi, phic):
        """Critical porosity model according to Nur’s modified Voigt average.

        Parameters
        ----------
        K0 : float or array-like
            mineral bulk modulus in GPa
        G0 : float or array-like
            mineral shear modulus in GPa
        phi : float or array-like
            porosity in frac
        phic : float
            critical porosity in frac

        Returns
        -------
        float or array-like
            K_dry,G_dry (GPa): dry elastic moduli of the framework
        """        
        
        K_dry = K0 * (1-phi/phic)
        G_dry = G0 * (1-phi/phic)

        return K_dry, G_dry
    
    @staticmethod
    def cripor_reuss(M0, Mf, phic, den=False):
        """In the suspension domain, the effective bulk and shear moduli of the rock can be estimated by using the Reuss (isostress) average.

        Parameters
        ----------
        M0 : float or array-like
            The solid phase modulus or density
        Mf : float or array-like
            The pore filled phase modulus or density
        phic : float
            critical porosity
        den : bool, optional
            If False: compute the reuss average for effective modulus of two mixing phases. If true, compute avearge density using mass balance, which corresponds to voigt average. Defaults to False.

        Returns
        -------
        float or array-like
            M (GPa/g.cc): average modulus or average density
        
        References
        ----------
        - Section 7.1 Rock physics handbook 2nd edition
        """        

        if den is False:

            M = EM.VRH(np.array([M0,Mf]), np.array([(1-phic,phic)]))[1]
        else:
            M = EM.VRH(np.array([M0,Mf]), np.array([(1-phic,phic)]))[0]

        return M

    @staticmethod
    def HS(f, K1, K2,G1, G2, bound='upper'):
        """Compute effective moduli of two-phase composite using hashin-strikmann bounds.

        Parameters
        ----------
        f : float
            0-1, volume fraction of stiff material
        K1 : float or array-like
            bulk modulus of stiff phase
        K2 : float or array-like
            bulk modulus of soft phase
        G1 : float or array-like
            shear modulus of stiff phase
        G2 : float or array-like
            shear modulus of soft phase
        bound : str, optional
            upper bound or lower bound. Defaults to 'upper'.

        Returns
        -------
        float or array-like
            K, G (GPa): effective moduli of two-phase composite
        """        

        if bound == 'upper':
            K=K1+ (1-f)/( (K2-K1)**-1 + f*(K1+4*G1/3)**-1 )

            Temp = (K1+2*G1)/(5*G1 *(K1+4*G1/3))
            G=G1+(1-f)/( (G2-G1)**-1 + 2*f*Temp)
        else:
            K=K2+ f/( (K1-K2)**-1 + (1-f)*(K2+4*G2/3)**-1 )

            Temp = (K2+2*G2)/(5*G2 *(K2+4*G2/3))
            G=G2+f/( (G1-G2)**-1 + 2*(1-f)*Temp)
        return K, G

    @staticmethod
    def Eshelby_Cheng(K, G, phi, alpha, Kf, mat=False):
        """Compute the effective anisotropic moduli of a cracked isotropic rock with single set fracture using Eshelby–Cheng Model.

        Parameters
        ----------
        K : float
            bulk modulus of the isotropic matrix GPa
        G : float
            shear modulus of the isotropic matrix GPa
        phi : float
            (crack) porosity
        alpha : float
            aspect ratio of crack
        Kf : float
            bulk modulus of the fluid. For dry cracks use fluid bulk modulus 0
        mat : bool, optional
            If true: the output is in matrix form, otherwise  is numpy array. Defaults to False.

        Returns
        -------
        _type_
            C_eff: effective moduli of cracked, transversely isotropic rocks

        References
        ----------
        - section 4.14 in The Rock Physics Handbook
        """        

        lamda = K-2*G/3

        sigma = (3*K-2*G)/(6*K+2*G)
        R = (1-2*sigma)/(8*np.pi*(1-sigma))
        Q = 3*R/(1-2*sigma)
        Sa = np.sqrt(1-alpha**2)
        Ia = 2*np.pi*alpha*(np.arccos(alpha)-alpha*Sa)/Sa**3
        Ic = 4*np.pi-2*Ia
        Iac = (Ic-Ia)/(3*Sa**2)
        Iaa = np.pi-3*Iac/4
        Iab = Iaa/3

        S11 = Q*Iab+R*Ia
        S33 = Q*(4*np.pi/3 - 2*Iac*alpha**2)+Ic*R
        S12 = Q*Iab-R*Ia
        S13 = Q*Iac*alpha**2-R*Ia
        S31 = Q*Iac-R*Ic
        S1212 = Q*Iab+R*Ia
        S1313 = Q*(1+alpha**2)*Iac/2 + R*(Ia+Ic)/2

        C = Kf/( 3*(K-Kf))
        D = S33*S11+S33*S12-2*S31*S13-(S11+S12+S33-1-3*C) - C*(S11+S12+2*(S33-S13-S31))
        E = S33*S11 - S31*S13 - (S33+S11-2*C-1) + C*(S31+S13-S11-S33)

        C11 = lamda*(S31-S33+1) + 2*G*E/ (D*(S12-S11+1))
        C33 = ((lamda+2*G)*(-S12-S11+1)+2*lamda*S13+4*G*C)/D
        C13 = ((lamda+2*G)*(S13+S31)-4*G*C+lamda*(S13-S12-S11-S33+2))/(2*D)
        C44 = G/(1-2*S1313)
        C66 = G/(1-2*S1212)

        if mat==False:
            C_eff = np.array([C11,C33,C13,C44,C66])
        else:
            C_eff = utils.write_VTI_matrix(C11,C33,C13,C44,C66)
        return C_eff

    # def Backus(V,lamda, G ):
    #     """Compute stiffnesses of a layered medium composed of thin isotropic layers using backus average model.

    #     Parameters
    #     ----------
    #     V : float or array-like
    #         volumetric fractions of N isotropic layering materials
    #     lamda : float or array-like
    #         Lamé coefficients of N isotropic layering materials
    #     G : float or array-like
    #         shear moduli of N isotropic layering materials

    #     Returns
    #     -------
    #     float or array-like
    #         C11,C33,C13,C44,C66: Elastic moduli of the anisotropic layered media
    #     """        

    #     C33=np.dot(V, 1/(lamda+2*G)) **-1
    #     C44=np.dot(V, 1/G)**-1
    #     C66=np.dot(V, G)
    #     C13=np.dot(V, 1/(lamda+2*G)) **-1 * np.dot(V, lamda/(lamda+2*G))
    #     C11=np.dot(V, 4*G*(lamda+G)/(lamda+2*G))+np.dot(V, 1/(lamda+2*G))**-1 * np.dot(V, lamda/(lamda+2*G))**2

    #     return C11,C33,C13,C44,C66


    @staticmethod
    def hudson(K, G, Ki, Gi, alpha, crd, order=1, axis=3):
        """Hudson’s effective crack model assuming weak inclusion for media with single crack set with all normals aligned along 1 or 3-axis. First and Second order corrections are both implemented. Notice that the second order correction has limitation. See Cheng (1993).

        Parameters
        ----------
        K : float 
            bulk modulus of isotropic background
        G : float 
            shear modulus of isotropic background
        Ki : float 
            bulk modulus of the inclusion material. For dry cracks: Ki=0
        Gi : float 
            shear modulus of the inclusion material
        alpha : float 
            crack aspect ratio
        crd : float 
            crack density
        order : int, optional
            approximation order.
                1: Hudson's model with first order correction.
                2: Hudson's model with first order correction.
                Defaults to 1.
        axis : int, optional
            axis of symmetry.
                1: crack normals aligned along 1-axis, output HTI
                3: crack normals aligned along 3-axis, output VTI
                Defaults to 3

        Returns
        -------
        _type_
            C_eff: effective moduli in 6x6 matrix form.
        """        
  
        lamda = K-2*G/3
        kapa = (Ki+4*Gi/3)*(lamda+2*G)/(np.pi*alpha*G*(lamda+G))
        M = 4*Gi*(lamda+G)/(np.pi*alpha*G*(3*lamda+4*G))
        U1 = 16*(lamda+2*G)/(3*(3*lamda+4*G)*(1+M))
        U3 = 4*(lamda+2*G)/(3*(lamda+G)*(1+kapa))
        # first order corrections
        C11_1 = -lamda**2*crd*U3/G
        C13_1 = -lamda*(lamda+2*G)*crd*U3/G
        C33_1 = -(lamda+2*G)**2*crd*U3/G
        C44_1 = -G*crd*U1
        #C66_1 = 0
        # second order corrections
        q = 15*lamda**2/G**2+28*lamda/G + 28
        C11_2 = q/15 * lamda**2/(lamda+2*G) *(crd*U3)**2
        C13_2 = q/15 * lamda*(crd*U3)**2
        C33_2 = q/15 * (lamda+2*G)*(crd*U3)**2
        C44_2 = 2/15 * G*(3*lamda+8*G)/(lamda+2*G)*(crd*U1)**2
        #C66_2 = 0
        if order == 1:
            C11 = lamda+2*G+C11_1
            C13 = lamda    +C13_1
            C33 = lamda+2*G+C33_1
            C44 = G        +C44_1
            C66 = G
        elif order == 2:
            C11 = lamda+2*G+C11_1+C11_2
            C13 = lamda    +C13_1+C13_2
            C33 = lamda+2*G+C33_1+C33_2
            C44 = G        +C44_1+C44_2
            C66 = G

        if axis ==3: # VTI
            C_eff = utils.write_VTI_matrix(C11,C33,C13,C44,C66)
        elif axis ==1: # HTI
            C_eff = utils.write_HTI_matrix(C33,C11,C13,C66,C44)
        return C_eff

    @staticmethod
    def hudson_rand(K, G, Ki, Gi, alpha, crd):
        """Hudson's crack model of a material containing randomly oriented inclusions. The model results agree with the consistent results of Budiansky and O’Connell (1976).

        Parameters
        ----------
        K : float or array-like
            bulk modulus of isotropic background
        G : float or array-like
            shear modulus of isotropic background
        Ki : float 
            bulk modulus of the inclusion material. For dry cracks: Ki=0
        Gi : float 
            shear modulus of the inclusion material, for fluid, Gi=0
        alpha : float 
            crack aspect ratio
        crd : float 
            crack density

        Returns
        -------
        float or array-like
            K_eff, G_eff (GPa): effective moduli of the medium with randomly oriented inclusions
        """        
        
        lamda = K-2*G/3
        kapa = (Ki+4*Gi/3)*(lamda+2*G)/(np.pi*alpha*G*(lamda+G))
        M = 4*Gi*(lamda+G)/(np.pi*alpha*G*(3*lamda+4*G))
        U1 = 16*(lamda+2*G)/(3*(3*lamda+4*G)*(1+M))
        U3 = 4*(lamda+2*G)/(3*(lamda+G)*(1+kapa))
        # first order corrections
        G_1 = -2*G*crd*(3*U1+2*U3)
        lamda_1 = (1/3)*(-(3*lamda+2*G)**2*crd*U3/(3*G) - 2*G_1)
        lamda_eff = lamda+lamda_1
        G_eff = G+G_1
        K_eff = lamda_eff+2*G_eff/3
        return K_eff,G_eff

    @staticmethod
    def hudson_ortho(K, G, Ki, Gi, alpha, crd):
        """Hudson’s first order effective crack model assuming weak inclusion for media with three crack sets with normals aligned along 1 2, and 3-axis respectively.  Model is valid for small crack density and aspect ratios.

        Parameters
        ----------
        K : float 
            bulk modulus of isotropic background
        G : float 
            shear modulus of isotropic background
        Ki : float 
            bulk modulus of the inclusion material. For dry cracks: Ki=0
        Gi : float 
            shear modulus of the inclusion material, for fluid, Gi=0
        alpha : nd array with size 3
            [alpha1, alpha2,alpha3] aspect ratios of  three crack sets
        crd : nd array with size 3
            [crd1, crd2, crd3] crack densities of three crack sets

        Returns
        -------
        2d array
            C_eff: effective moduli in 6x6 matrix form.
        """        
        
        # if type(alpha) == 'list' or type(crd) == 'list' :
        #     raise Exception("need python array as input for alpha and crd")
        alpha= np.asanyarray(alpha)
        crd= np.asanyarray(crd)
        lamda = K-2*G/3
        kapa = (Ki+4*Gi/3)*(lamda+2*G)/(np.pi*alpha*G*(lamda+G))
        M = 4*Gi*(lamda+G)/(np.pi*alpha*G*(3*lamda+4*G))
        U1 = 16*(lamda+2*G)/(3*(3*lamda+4*G)*(1+M))
        U3 = 4*(lamda+2*G)/(3*(lamda+G)*(1+kapa))

        # first order corrections
        C11_1 = -lamda**2*crd*U3/G
        C13_1 = -lamda*(lamda+2*G)*crd*U3/G
        C33_1 = -(lamda+2*G)**2*crd*U3/G
        C44_1 = -G*crd*U1
        C12_1 = C11_1 # C12= C11-2C66
        #C66_1 = 0
        C11=lamda+2*G+C33_1[0]+C11_1[1]+C11_1[2]
        C12=lamda    +C13_1[0]+C13_1[1]+C12_1[2]
        C13=lamda    +C13_1[0]+C12_1[1]+C13_1[2]
        C22=lamda+2*G+C11_1[0]+C33_1[1]+C11_1[2]
        C23=lamda    +C12_1[0]+C13_1[1]+C13_1[2]
        C33=lamda+2*G+C11_1[0]+C11_1[1]+C33_1[2]
        C44=G        +C44_1[1]+C44_1[2]
        C55=G        +C44_1[0]+C44_1[2]
        C66=G        +C44_1[0]+C44_1[1]

        C_eff = utils.write_matrix(C11,C22,C33,C12,C13,C23,C44,C55,C66)
        return C_eff

    @staticmethod
    def hudson_cone(K, G, Ki, Gi, alpha, crd, theta):
        """Hudson’s first order effective crack model assuming weak inclusion for media with crack normals randomly distributed at a fixed angle from the TI symmetry axis 3 forming a cone;

        Parameters
        ----------
        K : float 
            bulk modulus of isotropic background
        G : float 
            shear modulus of isotropic background
        Ki : float 
            bulk modulus of the inclusion material. For dry cracks: Ki=0
        Gi : float 
            shear modulus of the inclusion material, for fluid, Gi=0
        alpha : float 
            aspect ratios of crack sets
        crd : float 
            total crack density
        theta : float 
            the fixed angle between the crack normam and the symmetry axis x3. degree unit.

        Returns
        -------
        2d array
            C_eff: effective moduli of TI medium in 6x6 matrix form.
        """        
        
        theta= np.deg2rad(theta)
        lamda = K-2*G/3
        kapa = (Ki+4*Gi/3)*(lamda+2*G)/(np.pi*alpha*G*(lamda+G))
        M = 4*Gi*(lamda+G)/(np.pi*alpha*G*(3*lamda+4*G))
        U1 = 16*(lamda+2*G)/(3*(3*lamda+4*G)*(1+M))
        U3 = 4*(lamda+2*G)/(3*(lamda+G)*(1+kapa))

        # first order corrections
        C11_1 = -crd/(2*G)*(U3*(2*lamda**2+4*lamda*G*np.sin(theta)**2+3*G**2*np.sin(theta)**4)+U1*G**2*np.sin(theta)**2*(4-3*np.sin(theta)**2))

        C33_1 = -crd/G*(U3*(lamda+2*G*np.cos(theta)**2)**2+U1*G**2*4*np.cos(theta)**2*np.sin(theta)**2)

        C12_1 = -crd/(2*G)*(U3*(2*lamda**2+4*lamda*G*np.sin(theta)**2+G**2*np.sin(theta)**4)-U1*G**2*np.sin(theta)**4)

        C13_1 = -crd/G*(U3*(lamda+G*np.sin(theta)**2)*(lamda+2*G*np.cos(theta)**2)-U1*G**2*2*np.sin(theta)**2*np.cos(theta)**2)

        C44_1 = -crd/2*G*(U3*4*np.sin(theta)**2*np.cos(theta)**2+U1*(np.sin(theta)**2+2*np.cos(theta)**2-4*np.sin(theta)**2*np.cos(theta)**2))

        C66_1 = -crd/2*G*(U3*4*np.sin(theta)**4+U1*np.sin(theta)**2*(2-np.sin(theta)**2))
        #C22_1 = C11_1
        #C23_1 = C13_1
        #C55_1 = C44_1
        C11=lamda+2*G+C11_1
        C12=lamda    +C12_1
        C13=lamda    +C13_1
        #C22=C11
        #C23=C13
        C33=lamda+2*G+C33_1
        C44=G        +C44_1
        #C55=C44
        C66=G        +C66_1

        C_eff = utils.write_matrix(C11,C11,C33,C12,C13,C13,C44,C44,C66)
        return C_eff

    @staticmethod
    def Berryman_sc(K,G,X,Alpha):
        """Effective elastic moduli for multi-component composite using Berryman's Consistent (Coherent Potential Approximation) method.See also: PQ_vectorize, Berryman_func

        Parameters
        ----------
        K : array-like
            1d array of bulk moduli of N constituent phases, [K1,K2,...Kn]
        G : array-like
            1d array of shear moduli of N constituent phases, [G1,G2,...Gn]
        X : array-like
            1d array of volume fractions of N constituent phases, [x1,...xn], Sum(X) = 1.
        Alpha : array-like
            aspect ratios of N constituent phases. Note that α <1 for oblate spheroids and α > 1 for prolate spheroids, α = 1 for spherical pores,[α1,α2...αn]

        Returns
        -------
        array-like
            K_sc,G_sc: Effective bulk and shear moduli of the composite
        """        
        K= np.asanyarray(K)
        G= np.asanyarray(G)
        X=np.asanyarray(X)
        Alpha=np.asanyarray(Alpha)
        K_sc,G_sc=  fsolve(EM.Berryman_func, (K.mean(), G.mean()), args = (K,G,X,Alpha))
        return K_sc,G_sc

    @staticmethod
    def PQ_vectorize(Km,Gm, Ki,Gi, alpha):
        """compute geometric strain concentration factors P and Q for prolate and oblate spheroids according to Berymann (1980).See also: Berryman_sc, Berryman_func

        Parameters
        ----------
        Km : float
            Shear modulus of matrix phase. For Berryman SC       approach, this corresponds to the effective moduli of the composite.
        Gm : float
            Bulk modulus of matrix phase. For Berryman SC approach, this corresponds to the effective moduli of the composite.
        Ki : array-like
            1d array of bulk moduli of N constituent phases, [K1,K2,...Kn]
        Gi : array-like
            1d array of shear moduli of N constituent phases, [G1,G2,...Gn]
        alpha : array-like
            aspect ratios of N constituent phases. Note that α <1 for oblate spheroids and α > 1 for prolate spheroids, α = 1 for spherical pores,[α1,α2...αn]

        Returns
        -------
        array-like
            P,Q (array): geometric strain concentration factors, [P1,,,Pn],[Q1,,,Qn]
        """        
        dim = Ki.size
        P = np.empty(dim)
        Q = np.empty(dim)
        theta = np.empty(dim)
        alpha_ = alpha # copy and modfiy
        alpha_[alpha==1]=0.999 # when alpha_ is 1, the f is nan

        theta[alpha_<1]=alpha_[alpha_<1]/(1.0 - alpha_[alpha_<1]**2)**(3.0/2.0) * (np.arccos(alpha_[alpha_<1]) - alpha_[alpha_<1]*np.sqrt(1.0 - alpha_[alpha_<1]**2))

        theta[alpha_>1]=alpha_[alpha_>1]/(alpha_[alpha_>1]**2-1)**(3.0/2.0) * ( alpha_[alpha_>1]*(alpha_[alpha_>1]**2-1)**0.5 -np.cosh(alpha_[alpha_>1])**-1)

        f= alpha_**2*(3.0*theta - 2.0)/(1.0 - alpha_**2)
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
        # find and replace P and Q with alpha==1
        P[alpha==1]=(Km+4*Gm/3)/(Ki[alpha==1]+4*Gm/3)
        kesai= Gm/6 *(9*Km+8*Gm)/(Km+2*Gm)
        Q[alpha==1]= (Gm+kesai)/(Gi[alpha==1]+kesai)
        return P, Q

    @staticmethod
    def Berryman_func(params, K,G,X,Alpha ):
        """Form the system of equastions to solve. See 4.11.14 and 4.11.15 in Rock physics handbook 2020. See also: Berryman_sc

        Parameters
        ----------
        params : 
            Parameters to solve, K_sc, G_sc
        K : array
            1d array of bulk moduli of N constituent phases, [K1,K2,...Kn]
        G : array
            1d array of shear moduli of N constituent phases, [G1,G2,...Gn]
        X : array
            1d array of volume fractions of N constituent phases, [x1,...xn]
        Alpha : array
            aspect ratios of N constituent phases. Note that α <1 for oblate spheroids and α > 1 for prolate spheroids, α = 1 for spherical pores,[α1,α2...αn]

        Returns
        -------
        equation
            Eqs to be solved
        """        
        
        K_sc, G_sc = params
        P, Q = EM.PQ_vectorize(K_sc,G_sc, K,G, Alpha)
        eq1 = np.sum(X*(K-K_sc)*P)
        eq2 = np.sum(X*(G-G_sc)*Q)
        return  [eq1,eq2]


    @staticmethod
    def Swiss_cheese(Ks,Gs,phi): # Dilute_incl
        """Compute effective elastic moduli via "Swiss cheese" model with spherical pores. "Swiss cheese" model assumes a dilute distribution of spherical inclusions embedded in an * *unbounded* * homogenous solid.  It takes the "noninteracting assumption" in which all cavities (pores) are independent so that their contributions can be added.

        Parameters
        ----------
        Ks : float 
            Bulk modulus of matrix in GPa
        Gs : float 
            Shear modulus of matrix in GPa
        phi : float or array-like
            porosity

        Returns
        -------
        float or array-like
            Kdry,Gdry (GPa): effective elastic moduli
        """    
        Kdry=(1/Ks *(1+(1+3*Ks/(4*Gs))*phi))**-1
        Gdry=(1/Gs * (1+(15*Ks+20*Gs)*phi/(9*Ks+8*Gs)))**-1
        return Kdry, Gdry

    @staticmethod
    def SC(phi,Ks,Gs,iter_n):
        """Self-Consistent(SC) model with spherical pores considering the critical porosity and the interaction effect between inclusions.

        Parameters
        ----------
        phi : float or array-like
            porosity in frac, note that phi.shape== Ks.shape
        Ks : float
            bulk modulus of matrix phase in GPa
        Gs : float
            shear modulus of matrix phase in GPa
        iter_n : int
            iterations, necessary iterations increases as f increases.

        Returns
        -------
        float or array-like
            K_eff,G_eff (GPa): effective elastic moduli
        """     
        K_eff=Ks
        G_eff=Gs
        for i in range(iter_n):
            K_eff = (1/Ks + (1/K_eff+3/(4*G_eff))*phi) **-1
            G_eff= (1/Gs + (15*K_eff+20*G_eff)/(9*K_eff+8*G_eff) *phi/G_eff )**-1
        return K_eff,G_eff

    @staticmethod
    def Dilute_crack(Ks,Gs,cd):
        """The non-iteracting randomly oriented crack model.

        Parameters
        ----------
        Ks : float
            bulk modulus of uncracked medium in GPa
        Gs : float
            shear modulus of uncracked medium in GPa
        cd : float or array-like
            crack density

        Returns
        -------
        float or array-like
            K_eff,G_eff (GPa): effective elastic moduli
        """    
    
        nu= (3*Ks - 2*Gs)/(2*(3*Ks + Gs))  
        K_eff = ( 1/Ks * (1+16/9 *(1-nu**2)/(1-2*nu)*cd) )**-1
        G_eff = ( 1/Gs * (1+32*(1-nu)*(5-nu)/ (45*(2-nu)) * cd) )**-1

        return K_eff, G_eff 

    @staticmethod
    def OConnell_Budiansky(K0,G0,crd):
        """O’Connell and Budiansky (1974) presented equations for effective bulk and shear moduli of a cracked medium with randomly oriented dry penny-shaped cracks (in the limiting case when the aspect ratio α goes to 0)

        Parameters
        ----------
        K0 : float
            bulk modulus of background medium
        G0 : float
            shear modulus of background medium
        crd : float
            crack density

        Returns
        -------
        float
            K_dry,G_dry: dry elastic moduli of cracked medium
        """        
    
        nu0 = (3*K0-2*G0)/(6*K0+2*G0)#  Poisson ratio of the uncracked solid

        nu_eff = nu0*(1-16*crd/9) # approximation of the effective poisson'ratio of cracked solid
        K_dry = K0*(1-16*(1-nu_eff**2)*crd/(9*(1-2*nu_eff)))
        G_dry = G0*(1-32*(1-nu_eff)*(5-nu_eff)*crd/(45*(2-nu_eff)))
        return K_dry,G_dry


    @staticmethod
    def OConnell_Budiansky_fl(K0,G0,Kfl,crd, alpha):
        """Saturated effective elastic moduli using the O’Connell and Budiansky Consistent (SC) formulations under the constraints of small aspect ratio cracks with soft-fluid saturation.

        Parameters
        ----------
        K0 : float
            bulk modulus of background medium
        G0 : float
            shear modulus of background medium
        Kfl : float
            bulk modulus of soft fluid inclusion, e.g gas
        crd : float
            crack density
        alpha : float
            aspect ratio

        Returns
        -------
        float
            K_sat,G_sat: elastic moduli of cracked background fully saturated by soft fluid.

        References
	    ----------
        - O’Connell and Budiansky, (1974)
        """        
        
        nu0 = (3*K0-2*G0)/(6*K0+2*G0)#  Poisson ratio of the uncracked solid
        w = Kfl/alpha/K0
        # given crack density and w, solve for the D and nu_eff simulaneously using equations 23 and 25 in O’Connell and Budiansky, (1974)
        nu_eff, D =  fsolve(EM.OC_R_funcs, (0.2, 0.9), args = (crd,nu0,w))

        K_sat = K0*(1-16*(1-nu_eff**2)*crd*D/(9*(1-2*nu_eff)))
        G_sat = G0*(1-32/45*(1-nu_eff) *(D+ 3/(2-nu_eff))*crd)
        return K_sat,G_sat

    @staticmethod
    def OC_R_funcs(params, crd,nu_0,w ): # crd, nu_0,w
        """Form the system of equastions to solve. Given crack density and w, solve for the D and nu_eff simulaneously using equations 23 and 25 in O’Connell and Budiansky, (1974)

        Parameters
        ----------
        params : 
            Parameters to solve
        crd : float
            crack density
        nu_0 : float
            Poisson's ratio of background medium
        w : float
            softness indicator of fluid filled crack, w=Kfl/alpha/K0, soft fluid saturation is w is the order of 1

        Returns
        -------
        equation
            eqs to be solved
        """        
    
        nu_eff, D = params

        eq1 = 45/16 * (nu_0-nu_eff)/(1-nu_eff**2)*(2-nu_eff)/(D*(1+3*nu_0)*(2-nu_eff)-2*(1-2*nu_0)) - crd   # eq 23 in OC&R, 1974
        eq2 = crd * D**2-( crd+9/16 * (1-2*nu_eff)/(1-nu_0**2) + 3*w/(4*np.pi))*D + 9/16 * (1-2*nu_eff)/(1-nu_0**2) # eq 23 in OC&R, 1974
        return  [eq1,eq2]

    @staticmethod
    def PQ(Km,Gm, Ki,Gi, alpha): 
        """compute geometric strain concentration factors P and Q for prolate and oblate spheroids according to Berymann (1980). See also PQ_vectorize

        Parameters
        ----------
        Km : float
            Bulk modulus of matrix phase
        Gm : float
            Shear modulus of matrix phase
        Ki : float
            Bulk modulus of inclusion phase
        Gi : float
            Shear modulus of inclusion phase
        alpha : float
            aspect ratio of the inclusion. Note that α <1 for oblate spheroids and α > 1 for prolate spheroids

        Returns
        -------
        float
            P,Q (unitless): geometric strain concentration factors
        """       
        
        if alpha==1:
            P= (Km+4*Gm/3)/(Ki+4*Gm/3)
            kesai= Gm/6 *(9*Km+8*Gm)/(Km+2*Gm)
            Q= (Gm+kesai)/(Gi+kesai)
            
        
        else:

            if alpha<1:
                theta= alpha/(1.0 - alpha**2)**(3.0/2.0) * (np.arccos(alpha) - alpha*np.sqrt(1.0 - alpha**2))
            else:
                theta= alpha/(alpha**2-1)**(3.0/2.0) * ( alpha*(alpha**2-1)**0.5 -np.cosh(alpha)**-1)
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

    @staticmethod
    def DEM(y,t, params):
        '''
        ODE solver tutorial: https://physics.nyu.edu/pine/pymanual/html/chap9/chap9_scipy.html. 
        '''
        K_eff,G_eff=y  # unpack current values of y
        Gi,Ki,alpha = params # unpack parameters 
        P, Q= EM.PQ(G_eff,K_eff,Gi,Ki, alpha)
        derivs = [1/(1-t) * (Ki-K_eff) * P,  1/(1-t) * -G_eff * Q]
        return derivs

    @staticmethod
    def Berryman_DEM(Km,Gm, Ki, Gi, alpha,phi):
        """Compute elastic moduli of two-phase composites by incrementally adding inclusions of one phase (phase 2) to the matrix phase using Berryman DEM theory

        Parameters
        ----------
        Km : float
            host mineral bulk modulus 
        Gm : float
            host mineral shear modulus 
        Ki : float
            bulk modulus of inclusion
        Gi : float
            shear modulus of inclusion
        alpha : float
            aspect ratio of the inclusion phase
        phi : float
            desired fraction occupied by the inclusion
        """    
        #Bundle parameters for ODE solver
        params = [Gi,Ki,alpha]
        #Bundle initial conditions for ODE solver
        y0 = [Km,Gm]
        # Make time array for solution
        tStop = phi
        tInc = 0.01
        t = np.arange(0, tStop+tInc, tInc)
        psoln = odeint(EM.DEM, y0, t, args=(params,))
        K_dry_dem=psoln[:,0]
        G_dry_dem=psoln[:,1]
        return K_dry_dem, G_dry_dem,t

    @staticmethod
    def SC_dilute(Km, Gm, Ki, Gi, f, mode):
        """Elastic solids with elastic micro-inclusions. Random distribution of dilute spherical micro-inclusions in a two phase composite.

        Parameters
        ----------
        Km : float
            bulk modulus of matrix
        Gm : float
            shear modulus of matrix 
        Ki : float
            bulk modulus of inclusion 
        Gi : float
            shear modulus of inclusion
        f : float or array
            _descripvolume fraction of inclusion phase tion_
        mode : string
            'stress' if macro stress is prescribed. 'strain' if macro strain is prescribed.

        References
        ----------
        S. Nemat-Nasser and M. Hori (book) : Micromechanics: Overall Properties of Heterogeneous Materials. Sec 8

        Returns
        -------
        float or array
            K, G: effective moduli of the composite 
        """    
        
        nu = 0.5* (3*Km- 2*Gm)/(3*Km+Gm) # Poisson's ratio
        S1 = 1/3 * (1+nu)/( 1-nu)
        S2 = 2/15 * (4-5*nu)/(1-nu)
        if mode == 'stress':
            K = ( 1+ f * (Km/ (Km-Ki)- S1)**-1)**-1
            G = ( 1+ f * (Gm/ (Gm-Gi)- S2)**-1)**-1
        elif mode == 'strain':
            K = 1-f*(Km/(Km-Ki)-S1)**-1
            G = 1-f*(Gm/(Gm-Gi)-S2)**-1
        return K, G

    @staticmethod
    def SC_flex(f, iter_n, Km, Ki, Gm, Gi):
        """iteratively solving self consistent model for a two phase compposite consisting random distribution of spherical inclusion, not limited to pore. 

        Parameters
        ----------
        f : float or array
            volumetric fraction, f.shape== Km.shape
        iter_n : int
            iterations, necessary iterations increases as f increases.
        Km : float
            bulk modulus of matrix phase
        Ki : float
            bulk modulus of inclusion phase
        Gm : float
            shear modulus of matrix phase
        Gi : float
            shear modulus of inclusion phase 

        Reference
        ---------
        S. Nemat-Nasser and M. Hori (book) : Micromechanics: Overall Properties of Heterogeneous Materials. Sec 8

        Returns
        -------
        float or array
            K_eff, G_eff (GPa): Effective elastic moduli
        """    
        
        K_eff=Km
        G_eff=Km
        for i in range(iter_n):
            nu = 0.5* (3*K_eff- 2*G_eff)/(3*K_eff+G_eff) # Poisson's ratio
            S1 = 1/3 * (1+nu)/( 1-nu)
            S2 = 2/15 * (4-5*nu)/(1-nu)
            K_eff = (1- f* K_eff*(Km-Ki)/(Km*(K_eff-Ki)) *(K_eff/(K_eff-Ki)-S1)**-1) * Km

            G_eff= (1- f * G_eff*(Gm-Gi)/(Gm*(G_eff-Gi)) *(G_eff/(G_eff-Gi)-S2)**-1 ) *Gm 
        return K_eff,G_eff

    @staticmethod
    def MT_average(f, Kmat, Gmat, K1, G1, K2, G2): 
        """Compute Two-phase composite without matrix using modified Mori-Takana Scheme according to  Iwakuma 2003, one of the inhomogeneities must be considered as a matrix in the limiting model.

        Parameters
        ----------
        f : float or array
            Volume fraction of matrix/inhomogeneity 1. f1=1-f2, (1-f) can be regarded as pseudo crack density.
        Kmat : float
            Bulk modulus of matrix/inhomogeneity 1
        Gmat : float
            shear modulus of  matrix/ inhomogeneity 1
        K1 : float
            Bulk modulus of inhomogeneity 1
        G1 : float
            shear modulus of inhomogeneity 1
        K2 : float
            Bulk modulus of inhomogeneity 2
        G2 : float
            shear modulus of inhomogeneity 2

        Returns
        -------
        float or array
            K_ave, G_ave [GPa]: MT average bulk and shear modulus 
        """    
    
        possion= (3*Kmat -2*Gmat)/(2*(3*Kmat+Gmat)) # possion (unitless): possion's ratio of virtual matrix 
        
        alpha=(1+possion) / ( 3*(1-possion) )
        beta= 2* (4-5*possion) / ( 15* (1-possion) )

        K_ave= (f*Kmat*K1/(Kmat-(Kmat-K1)*alpha) + (1-f)*Kmat*K2/(Kmat-(Kmat-K2)*alpha)) / ( f*Kmat/(Kmat-(Kmat-K1)*alpha) + (1-f)*Kmat/(Kmat-(Kmat-K2)*alpha)  )
        G_ave= (f*Gmat*G1/(Gmat-(Gmat-G1)*beta) + (1-f)*Gmat*G2/(Gmat-(Gmat-G2)*beta)) / ( f*Gmat/(Gmat-(Gmat-G1)*beta) + (1-f)*Gmat/(Gmat-(Gmat-G2)*beta)  )
        return K_ave, G_ave