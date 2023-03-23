#!/usr/bin/env python
# -*-coding:utf-8 -*-


import numpy as np


class Empirical:
    """ Empirical relations that widely applied 
    """   
     
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