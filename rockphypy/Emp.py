#!/usr/bin/env python
# -*-coding:utf-8 -*-


import numpy as np


class Empirical:
    """ Empirical relations that widely applied 
    """   
    @staticmethod   
    def krief(phi, Kg,Gg):
        """Compute porous background elastic constants as a function of porosity according to Krief model. 

        Parameters
        ----------
        phi : float or array
            porosity in porous rock 
        Kg : float
            grain bulk modulus 
        Gg : float
            grain shear modulus 

        Returns
        -------
        float or array
            lamda,G,K
        """        
    
        a0= 1-(1-phi)**(3/(1-phi))
        G = (1-a0)*Gg
        K = (1-a0)*Kg
        lamda= K-2*G/3
        return lamda,G,K

    @staticmethod   
    def esti_VS(Vp, vsh):
        """Estimate, using the Greenberg-Castagna empirical relations, the shearwave velocity in a brine-saturated shaly sandstone with vp
        since we only assume two minearl phases: so L=2, X1= 1-vsh, X2= vsh

        Parameters
        ----------
        Vp : float or array-like
           compressional velocities, m/s
        vsh : float or array-like
            shale volume

        Returns
        -------
        float or array-like
            vs (m/s): estimated shear wave velocities 

        References
	    ----------
        - handbook of rock physics P516
        """           
        
        Vs_sand= 0.80416*Vp/1000 - 0.85588 # km/s
        Vs_shale= 0.76969*Vp/1000 - 0.86735 #km/s 

        Vs_arith= (1-vsh)*Vs_sand+ vsh*Vs_shale
        Vs_harm=( (1-vsh)/Vs_sand +vsh/Vs_shale )**-1
        Vs= 0.5*(Vs_arith+Vs_harm)
        return Vs*1000

    @staticmethod
    def han(phi,C):
        """Han (1986) found empirical regressions relating ultrasonic (laboratory) velocities to porosity and clay content.effective pressure is 20Mpa

        Parameters
        ----------
        phi : float or array-like
            porosity
        C : float or array-like
            clay volume fraction
        Returns
        -------
        float or array-like
            P and S wave velocities
        """        
             
        vp=5.59-6.93*phi-2.81*C
        vs=3.52-4.91*phi-1.89*C
        return vp*1000,vs*1000
########################################################## to do 
    @staticmethod
    def ehrenberg(Z):
        """porosity reference trend for Norwegian Sea sandstone. Note that the functional form of the porosity model is not published in Ehrenberg (1990). It is obtained by linear regression of the digitized data point from the original plot in the paper. 

        Parameters
        ----------
        Z : float or array
            burial depth below see floor in Km

        References
        ----------
        Ehrenberg, S., 1990, Relationship between diagenesis and reservoir quality in sandstones of the Garn Formation, Haltenbanken, mid-Norwegian continental shelf: AAPG bulletin, 74, no. 10, 1538

        Returns
        -------
        float or array
            porosity
        """        
           
        Z= Z/1000
        phi = -0.092200747*Z+0.4802422
        return phi 

    @staticmethod
    def yu_segment_trend(Z):
        """Reference trend for Norwegian sea normally buried clean sandstones 

        Parameters
        ----------
        Z : float or array
            burial depth below see floor in m

        Returns
        -------
        float or array
            P wave velocities
        """        
        
        V1=1707.998+0.66*Z[Z<=2630]
        V2=1200.84+0.853*Z[Z>2630]
        V= np.hstack((V1,V2))
        return V

    @staticmethod
    def ramm_porosity(Z, HB=True):
        """porosity reference trend according to Ramm & Bjørlykke (1994) 

        Parameters
        ----------
        Z : float or array
            burial depth wrt. sea floor in m
        HB : bool, optional
            if True: only show the regression line for halten bakken area porosity data False: The regression line for all porosity from north sea and norwegian sea, by default True

        Returns
        -------
        float or array
            porosity
        """        
        
        if HB==True:
            phi=  46.4-0.0085*Z
        else:
            phi= 42.7-0.0069*Z
        return phi/100

    @staticmethod
    def ramm_porosity_segment(Z):
        """segment porosity reference trend according to Ramm & Bjørlykke (1994) considering the mechanical and chemical compaction

        Parameters
        ----------
        Z : float or array
            burial depth wrt. sea floor in m

        Returns
        -------
        float or array
            porosity
        """        
        phi1= 45*np.exp(-0.00025*Z)
        phi2= 25-0.013*(Z-2500)
        phi1[phi1>phi2]=phi2[phi2<phi1]
        return phi1

    @staticmethod
    def empirical_StPeter(Pe, sample=1):
        """compute the Vp and Vs for st peter sandstone using the empirical relationship in the form of V= A+KPe-Be^(-DPe)

        Parameters
        ----------
        Pe : float or array
            effective pressure in Kbar, 1kbar= 100Mpa
        sample : int, optional
            1-sample 1, phi= 0.205. 2-sample 2, phi= 0.187, by default 1

        Returns
        -------
        float or array
            Vp,Vs in km
        """        
        
        if sample==1:

            Vp= 4.21+0.187*Pe-0.746*np.exp(-24*Pe)
            Vs= 2.58+0.160*Pe-0.781*np.exp(-20*Pe)

        elif sample ==2:
            Vp= 4.55+0.298*Pe-0.8*np.exp(-19*Pe)
            Vs= 2.83+0.195*Pe-0.741*np.exp(-16*Pe)
        return Vp*1000, Vs*1000

    @staticmethod
    def Scherbaum(Z):
        """ velocity depth trend for Lower and Middle Buntsandstein

        Parameters
        ----------
        Z : float or array
            burial depth wrt. sea floor in m

        References
        ----------
        Scherbaum, F., 1982. Seismic velocities in sedimentary rocks—indicators of subsidence and uplift?. Geologische Rundschau, 71(2), pp.519-536.

        Returns
        -------
        float or array
            P wave velocities in m/s
        """        
        V = 2325+0.51*Z
        return V

    @staticmethod
    def Sclater(phi):
        """Sclater-Christie exponential curve for sandstone

        Parameters
        ----------
        phi : float or array
            porosity

        Returns
        -------
        Z : float or array
            depth wrt. sea floor in km.
        """        
        Z= 3.7*np.log(0.49/phi)
        return Z

    @staticmethod
    def Storvoll(Z):
        """Storvoll velocity compaction trend. The trend is for shale and shaly sediments but also used for siliciclastic rock like sandstone 

        Parameters
        ----------
        Z : float or array
            depth wrt. sea floor in m.

        Returns
        -------
        float or array
            Vp meters per second
        """        
        V= 1/1.76*Z+2600
        return V

    @staticmethod
    def Hillis(Z):
        """compaction trend for Bunter Sandstone in North sea 

        Parameters
        ----------
        Z : float or array
            depth below sea bed (in kilometers)

        Returns
        -------
        float or array
            Vp km/s
        """        
        transit= 135.9-20.22*Z/1000
        V= 304.8*1000/transit
        return V

    @staticmethod
    def Japsen(Z):
        """a segmented linear velocity–depth function, These equations are considered as approximation for bunter sandstone trend although they are originally for bunter shale. proposed by Japsen 1999

        Parameters
        ----------
        Z : float or array
            depth below sea bed in m.

        Returns
        -------
        float or array
            Vp m/s
        """        
        
        V1=1550+0.6*Z[Z<1393]
        V2= -400+2*Z[(Z>=1393)&(Z<2000)]
        V3=2600+0.5*Z[(Z>=2000)&(Z<3500)]
        V4=3475+0.25*Z[(Z>=3500)&(Z<5300)]
        V= np.hstack((V1,V2,V3,V4))
        return V

    @staticmethod    
    def hjelstuen(Z):
        """Velocity-depth relationships for the Bjørna-Sørkapp margin deposits. note:  the seismic velocities are not directly comparable with velocities from sonic logs (because of the different frequencies), and the velocity-depth profile of Hjelstuen et al. (1996) has not been corrected for uplift and erosion

        Parameters
        ----------
        Z : float or array
            Z< 3.8km

        Returns
        -------
        float or array
            V: m/s
        """        
          
        V=1.87+0.55*(Z/1000) # convert Z to km 
        return V*1000

    @staticmethod
    def Cp(phi):
        """The coordination number n depends on porosity, as shown by
        Murphy 1982.

        Parameters
        ----------
        phi : float or array
            total porosity , for a porosity of 0.4, n=8.6

        Returns
        -------
        float or array
            coordination number
        """        
        return 20-34*phi+14*phi**2