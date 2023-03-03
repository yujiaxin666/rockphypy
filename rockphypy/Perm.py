#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   Perm.py
@Time    :   2023/02/19 22:11:45
@Author  :   Jiaxin Yu 
@Contact :   yujiaxin666@outlook.com


Recreation and modifcations of the Permeabilities models in Rock physics handbook matlab tools. 
'''

class Permeability:
    """_summary_
    """    
    def Kozeny_Carman(phi,d):
        """Describe the permeability in a porous medium using Kozeny-Carman equation assuming the turtuosity tau=sqrt(2), 1/B=2.5 for unconsolidated monomodal sphere pack. 

        Args:
            phi (frac): porosity
            d (m):pore diameter 
        Examples:
            phi= np.linspace(0.01,0.35,100)
            d= 250
            k= Kozeny_Carman(phi, d)
            plt.semilogy(phi, k )
        Returns:
            k (m^2): the resulting permeability is in the same units as d^2 
        """    
        k= d**2/180 *phi**3/(1-phi)**2

        return k

    def Kozeny_Carman_Percolation(phi,phic, d, B):
        """The Kozeny−Carman relations incorporating the percolation effect 

        Args:
            phi (frac): porosity
            phic (frac):percolation porosity
            d (m): pore diameter  
            B : geometric factor that partly accounts for the irregularities of pore shapes.

        Returns:
            k (m^2): the resulting permeability is in the same units as d^2 
        """    
        
        k= B*d**2* (phi-phic)**3/(1+phic-phi)**2

        return k

    def Owolabi(phi, Swi):
        """Estimate the permeability in uncosonlidated sands of Pleistocene to Oligocene age in Eastern Niger Delta from log derived porosityand irreducible water saturation.

        Args:
            phi (frac): porosity
            Swi (frac): irreducible water-saturation from welllogs


        Returns:
            k_oil, k_gas: permeabilities in mD for oil and gas sand reservoir, respectively
        """    
        

        #For Oil Sand
        k_oil = 307 + (26552*(phi**2))-(34540*(phi*Swi)**2)
        #For Gas Sand
        k_gas = 30.7 + (2655*(phi**2))-(3454*(phi*Swi)**2)
        return k_oil, k_gas 

    def Perm_logs(phi, Swi):
        """Various empirical correlations of between permeability, porosity and irreducible water-saturation from welllogs. Models includs Tixier, Timur, Coates and Coates-Dumanoir. 
        
        Assumption: 
        
            The functional forms used in these equations have to be calibrated, whenever possible, to site-specific data.
            The rock is isotropic.
            Fluid-bearing rock is completely saturated.

        Args:
            phi (frac): porosity
            Swi (frac): irreducible water-saturation from welllogs

        Returns:
            k_tixier, k_Timur , k_coates, k_coates_Dumanoir: different permeability estimations, in the unit of mD 
        """    

        k_tixier = 62500*phi**6/Swi**2

        k_Timur = 10000*phi**4.5/Swi**2

        k_coates = 10000*phi**4*(1-Swi)**2/Swi**2

        k_coates_Dumanoir = 352*phi**4/Swi**4
        
        return k_tixier, k_Timur , k_coates, k_coates_Dumanoir


    def Panda_Lake(d, C,S,tau,phi):
        """ Modified Kozeny-carman relation incorpating the contribution of grain size variation and sorting using Manmath N. Panda and Larry W. Lake relation. 

        Args:
            d (micro m):  mean particles size
            C (_type_): coefficient of variation of particles size distribution
            S (_type_): skewness of particles size distribution
            tau (_type_): tortuosity factor
            phi (frac): porosity
        Refs:
            Estimation of Single-Phase permeability from parameters of particle-Size Distribution, Manmath N. Panda and Larry W. Lake, AAPG 1994.
        Returns:
        k (md): permeability
        """    
        
        
        k= d**2*phi**3*(C**3*S+3*C**2+1)**2/(72*tau*(1-phi)**2*(C**2+1)**2)


        return k


    def Panda_Lake_cem(phi,d):
        """ Quantify the effects of cements on the single phase permeability estimate of unconsolidated sand using Panda & Lake model

        Args:
            phi (frac): porosity
            d (micro m):  mean particles size
        Returns:
        k (md): permeability
        """    
        K =3.34*(d**2)*((phi**3)/(1-phi)**2)
        return 


    def Revil(phi, d):
        """Estimate permeability in very shaly rock using Revil et al. 1997 

        Args:
            phi (frac): porosity
            d (micro m): average grain diameter

        Returns:
            k (md): permeability
        """    
        
        k = 1000*(d**2)*(phi**4.5)/24

        return k

    def Fredrich(phi, d, b):
        """ Compute permability considering Pore Geometry and Transport Properties of Fontainebleau Sandstone

        Args:
            phi (frac): porosity>10%
            d (micro m): pore diameter 
            b : shape factor b is equal to 2 for circular tubes and equal to 3 for cracks. 
        Refs:
            Fredrich, J. T., Greaves, K. H., & Martin, J. W. (1993, December). Pore geometry and transport properties of Fontainebleau sandstone. In International journal of rock mechanics and mining sciences & geomechanics abstracts (Vol. 30, No. 7, pp. 691-697). Pergamon.
        Returns:
            k (md): permeability
        """    

        # formation factor assuming tau^2=2.5
        F= 2.5/phi
        # pore surface area per unit sample volume
        Sv= 6*(1-phi)*d
        # permeability
        k= 1/(b*F) * (phi/Sv)**2
        return k

    def Bloch(S,C,D):
        """Predict porosity and permeability in sandstones prior to drilling using Bloch empirical relations obtain in Yacheng field.

        Args:
            S : Trask sorting coefficient
            C (frac): Rigid grain content
            D (mm): Grain size
        Returns:
            phi, k: porosity (frac) and permeability (mD), respectively 
        """   
        # Porosity Calculation
        phi = -6.1 + (9.8/S) + (0.17*C)
        # Permeability Calculation
        k = 10**(-4.67 + (1.34*D) + (4.08/S)+ 3.42*(C/100))
        return phi, k


    def Bernabe(phi, crf,w, r):
        """Bernabe models permit to compute the permeability and porosity of strongly pressure dependent pores such as cracks and approximately constant pores associated with tubes and nodal pores.
        Args:
            phi (frac): total porosity
            crf (frac): crack fraction in pore volume
            w (micro meter): width or aperture of the equivalent crack
            r (micro meter): radius of the tube 
        Refs:
            Bernabe, Y. (1991). Pore geometry and pressure dependence of the transport properties in sandstones. Geophysics, 56(4), 436-446.
        Returns:
            k (md): total permeability
        """    
        #Crack porosity
        phicrack = phi*crf
        #Crack permeability
        # Kcrack = (w^2)*Phicrack/(12*Taucrack^2); (w^2)*(crf*phi)/30 ,  Taucrack = crack tortuosity, Taucrack^2 = 2.5,     
        kcrack = (w**2)*phicrack/30
        #Tube porosity Ktube = (r^2)*Phitube/(8*Tautube^2); (r^2)*((1-CrackFraction)*Phi)/20; tube tortuosity Tautube^2 = 2.5
        phitube = phi - phicrack
        #Tube permeability 
        ktube = (r**2)*phitube/20
        #Total permeability
        k = kcrack+ktube
        return k 