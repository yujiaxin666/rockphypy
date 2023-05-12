#!/usr/bin/env python
# -*-coding:utf-8 -*-


import numpy as np

class Anisotropy:
    """Effective models, coordinate transform, anisotropic parameters and phase velocities that can be applied to anisotropic media.  
    """    

    @staticmethod
    def Thomsen(C11, C33, C13, C44, C66, den, theta):
        """Compute thomsen parameters and three phase velocities for weak anisotropic TI media with vertical symmetry axis. 

        Parameters
        ----------
        Cij : float
            Stiffnesses in GPa
        den : float
            density of the effective medium, for backus average, the effective density can be computed using VRH method for Voigt average. 
        theta : float or array-like
            angle of incidence

        Returns
        -------
        float or array-like
            VP, VSV, VSH: wave velocities propagating along given direction
        """        
            
        alpha= np.sqrt(C33/den) *1e3
        beta= np.sqrt(C44/den) *1e3
        epsilon= 0.5* (C11-C33)/C33
        gamma= 0.5 * (C66-C44)/C44
        theta= np.deg2rad(theta) # degree to radian
        delta= ((C13+C44)**2-(C33-C44)**2 )/ (2*C33*(C33-C44))
        VP= alpha*(1+delta*np.sin(theta)**2*np.cos(theta)**2+epsilon*np.sin(theta)**4)
        VSV= beta*(1+alpha**2/beta**2 * (epsilon-delta)*np.sin(theta)**2*np.cos(theta)**2)
        VSH= beta*(1+gamma*np.sin(theta)**2)
        return VP, VSV, VSH, epsilon,gamma,delta

    @staticmethod
    def Thomsen_Tsvankin(C11,C22,C33,C12,C13,C23,C44,C55,C66):
        """Elastic constants of an orthorhombic elastic medium defined by Tsvankin’s notation for weak elastic anisotropy assuming the vertical symmetry axis is along the x3 direction.

        Parameters
        ----------
        Cij : float
            Stiffnesses in GPa

        Returns
        -------
        floats
            Thomsen-Tsvankin parameters 
        """        
    
        # defined in the symmetry plane with normal in x1 direction);
        epsilon_1 = (C22-C33)/(2*C33) 

        delta_1 = ((C23+C44)**2-(C33-C44)**2)/(2*C33*(C33-C44))
        gamma_1 = (C66-C55)/2/C55
        # defined in the symmetry plane with normal in x2 direction);
        epsilon_2 = (C11-C33)/(2*C33) 
        delta_2 = ((C13+C55)**2-(C33-C55)**2)/(2*C33*(C33-C55))
        gamma_2 = (C66-C44)/2/C44

        #gamma_= (C44-C55)/(2*C55)

        # defined in the horizontal plane with normal in x3 direction
        delta_3 = ((C12+C66)**2-(C11-C66)**2)/2/C11/(C11-C66)

        return epsilon_1,delta_1, gamma_1, epsilon_2, delta_2,gamma_2, delta_3

    @staticmethod
    def Backus(V,lamda, G ):
        """Computes stiffnesses of a layered medium using backus average model. 

        Parameters
        ----------
        V : float or array-like
            volumetric fractions of N isotropic layering materials
        lamda : float or array-like
            Lamé coefficients of N isotropic layering materials
        G : float or array-like
            shear moduli of N isotropic layering materials

        Returns
        -------
        float or array-like
            C11,C33,C13,C44,C66:Elastic moduli of the anisotropic layered media
        """        
        V= np.asanyarray(V)
        V=V/np.sum(V) # normalize
        lamda= np.asanyarray(lamda)
        G= np.asanyarray(G)
        # if isinstance(V, ( np.ndarray)) is False:
        #     V=np.array(V)
        # V=V/np.sum(V) # normalize 
        # if isinstance(lamda, ( np.ndarray)) is False:
        #     lamda=np.array(lamda)
        # if isinstance(G, ( np.ndarray)) is False:
        #     G=np.array(G)
        
        C33=np.dot(V, 1/(lamda+2*G)) **-1
        C44=np.dot(V, 1/G)**-1
        C66=np.dot(V, G)
        C13=np.dot(V, 1/(lamda+2*G)) **-1 * np.dot(V, lamda/(lamda+2*G))
        C11=np.dot(V, 4*G*(lamda+G)/(lamda+2*G))+np.dot(V, 1/(lamda+2*G))**-1 * np.dot(V, lamda/(lamda+2*G))**2
        
        return C11,C33,C13,C44,C66

    @staticmethod
    def Backus_log(Vp,Vs,Den,Depth):
        """Computes Backus Average from log data, notice that the Depth is 1d Vector including each top depth of layer and also the bottom of last layer. 

        Parameters
        ----------
        Vp : array 
            P wave velocities of layers [Vp1,Vp2...Vpn], Km/s, size N
        Vs : array
            S wave velocities of layers [Vs1,Vs2...Vsn],Km/s size N
        Den : array
            Densities of layers, size N
        Depth : array
            1d depth, ATTENTION: each depth point corresponds to the top of thin isotropic layer, the bottom of the sedimentary package is the last depth point. [dep1,dep2,,,,,,depn, depn+1], size N+1

        Returns
        -------
        array-like
            Stiffness coeffs and averaged density
        """        
        Vp, Vs, Den, Depth= np.asanyarray(Vp),np.asanyarray(Vs), np.asanyarray(Den),np.asanyarray(Depth)
        # compute lame constants from log data 
        G = Den*Vs**2
        lamda = Den*Vp**2 - 2*G
        # compute thickness of each layer given depth  
        thickness = np.diff(Depth)
        # Volume fraction
        V= thickness/thickness.sum()
        
        C11,C33,C13,C44,C66 = Anisotropy.Backus(V,lamda, G )

        # Mass balance for density 
        den = np.dot(V, Den)

        return C11,C33,C13,C44,C66, den

    @staticmethod
    def vel_azi_HTI(C,Den,azimuth):
        """Given stiffnesses and density of the HTI medium, compute the azimuth dependent phase velocities.

        Parameters
        ----------
        C : 2d array
            stiffness matrix of the HTI medium
        Den : float
            density of the fractured medium
        azimuth : float or array like
            azimuth angle, degree

        Returns
        -------
        float or array like
            VP,VSH, VSV: phase velocities 
        """        
  
        azimuth= np.radians(azimuth)
        C11=C[0,0]
        C12=C[0,1]
        C33=C[2,2]
        C44=C[3,3]
        C55=C[4,4]
        R= np.sqrt( ((C33-C55)*np.sin(azimuth)**2- (C11-C55)*np.cos(azimuth)**2)**2 +4*(C12+C55)**2 *np.sin(azimuth)**2*np.cos(azimuth)**2  )
        Vp = 0.5/Den *( (C33+C55 )*np.sin(azimuth)**2 +(C11+C55)*np.cos(azimuth)**2 + R)
        Vsh= 0.5/Den *( (C33+C55 )*np.sin(azimuth)**2 +(C11+C55)*np.cos(azimuth)**2 - R)
        Vsv= 1/Den * (C55+ (C44-C55)*np.sin(azimuth)**2)
        VP= np.sqrt(Vp)
        VSH= np.sqrt(Vsh)
        VSV= np.sqrt(Vsv)
        return VP,VSH, VSV

    @staticmethod
    def vel_azi_VTI(C,Den,azimuth):
        """Given stiffnesses and density of the VTI medium, compute the azimuth dependent phase velocities.

        Parameters
        ----------
        C : 2d array
            stiffness matrix of the VTI medium
        Den : float
            density of the fractured medium
        azimuth : float or array like
            azimuth angle, degree

        Returns
        -------
        float or array like
            VP,VSH, VSV: phase velocities 
        """        
  
        azimuth= np.radians(azimuth)

        C11=C[0,0]
        C13=C[0,2]
        C33=C[2,2]
        C44=C[3,3]
        C66=C[5,5]

        R= np.sqrt( ((C11-C44)*np.sin(azimuth)**2- (C33-C44)*np.cos(azimuth)**2)**2 +(C13+C44)**2 *np.sin(2*azimuth)**2 )
        Vp = 0.5/Den *( C11*np.sin(azimuth)**2 +C33*np.cos(azimuth)**2 + C44 + R)
        Vsv= 0.5/Den *( C11*np.sin(azimuth)**2 +C33*np.cos(azimuth)**2 + C44 - R)
        Vsh= 1/Den * (C66*np.sin(azimuth)**2+C44*np.cos(azimuth)**2 )
        VP = np.sqrt(Vp)
        VSH= np.sqrt(Vsh)
        VSV= np.sqrt(Vsv)
        return VP,VSH, VSV

    @staticmethod
    def Bond_trans(C, theta, axis=3):
        """Coordinate Transformations for stiffness matrix in 6x6 Voigt notation using Bond transformation matrix. 

        Parameters
        ----------
        C : 2d array
            original stiffness matrix
        theta : float
            rotational angle
        axis : int, optional
            axis=1: fix 1-axis, rotate 2 and 3 axis, examples can be a TTI(Tilted TI) resulted from the rotation of VTI with horizontal aligned fracture sets wrt the vertical x3 axis. In this case, the input C should be a VTI matrix
            axis=3: fix 3-axis, rotate 1 and 2 axis, E.g. seismic measurements of HTI media e.g caused by vertically aligned fractures. The angle theta may be assigned to be the angle between the fracture normal and a seismic line.


        Returns
        -------
        2d array
            C_trans, new stiffness matrix wrt to the original right-hand rectangular Cartesian coordinates

        References
        ----------
        - Bond, W., Jan. 1943, The mathematics of the physical properties of crystals, The Bell System Technical Journal, 1-72.
        
        """        
 
        theta = np.deg2rad(theta)
        # axis 3 is fixed, rotate x1 and x2 axes.
        if axis==3:
            b11 = np.cos(theta) 
            b22 = np.cos(theta) 
            b33 = 1
            b21 = np.cos(theta + np.pi/2) 
            b12 = np.cos(np.pi/2 - theta)
            b13 = 0 
            b31 = 0
            b23 = 0 
            b32 = 0
        # axis 1 is fixed, rotate x2 and x3 axes.
        elif axis==1:
            b11 = 1
            b22 = np.cos(theta) 
            b33 = np.cos(theta) 
            b21 = 0
            b12 = 0
            b13 = 0 
            b31 = 0
            b23 = np.cos(theta + np.pi/2)  
            b32 = np.cos(np.pi/2 - theta)

        # Bond matrix M
        
        M1 = np.array([b11**2, b12**2, b13**2, b21**2, b22**2, b23**2, b31**2, b32**2, b33**2]).reshape((3,3))
        M2 = np.array([b12*b13, b13*b11, b11*b12,
            b22*b23, b23*b21, b21*b22, 
            b32*b33, b33*b31, b31*b32]).reshape((3,3))
        M3 = np.array([b21*b31, b22*b32, b23*b33,
            b31*b11, b32*b12, b33*b13,
            b11*b21, b12*b22, b13*b23]).reshape((3,3))
        M4 = np.array([b22*b33+b23*b32, b21*b33+b23*b31, b22*b31+b21*b32,
            b12*b33+b13*b32, b11*b33+b13*b31, b11*b32+b12*b31,
            b22*b13+b12*b23, b11*b23+b13*b21, b22*b11+b12*b21]).reshape((3,3))

        M = np.block([[M1, M2], [M3, M4]])
        C_trans= M @ C @ M.T # matrix multiplication
        return C_trans

