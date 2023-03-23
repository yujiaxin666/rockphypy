#!/usr/bin/env python
# -*-coding:utf-8 -*-


'''
Batzle and Wang functionalities

'''
import numpy as np
import warnings

class BW:
    """Effective CO2, natural gas, brine and oil property calculation using original and modified Batzle-Wang equations.
    """    
    def dz_dp(P_pr,T_pr):
        """Values for dZ/dPpr obtained from equation 10b in Batzle and Wang (1992).
        """    
        # analytic
        dzdp=(0.03+0.00527*(3.5-T_pr)**3)+0.109*(3.85-T_pr)**2*1.2*P_pr**0.2*-(0.45+8*(0.56-1/T_pr)**2)/T_pr*np.exp(-(0.45+8*(0.56-1/T_pr)**2)*P_pr**1.2/T_pr)

        # numerical approximation 
        # dzdp= 1.938783*P_pr**0.2*(1 - 0.25974025974026*T_pr)**2*(-8*(0.56 - 1/T_pr)**2 - 0.45)*np.exp(P_pr**1.2*(-8*(0.56 - 1/T_pr)**2 - 0.45)/T_pr)/T_pr + 0.22595125*(1 - 0.285714285714286*T_pr)**3 + 0.03
        return dzdp

    #-----------------Pure Gas--------------------------#
    def pseudo_p_t(P,T,G):
        """Calculate the pseudoreduced temperature and pressure according to Thomas et al. 1970.

        Parameters
        ----------
        P : float or array-like
            Pressure in MPa
        T : float or array-like
            Temperature in °C
        G : float or array-like
            Gas gravity

        Returns
        -------
        float or array-like
            Ta: absolute temperature
            Ppr:pseudoreduced pressure 
            Tpr:pseudoreduced temperature 
        """        
 
        # convert the temperature to absolute temperature
        Ta=T+273.15
        P_pr=P/(4.892-0.4048*G)
        T_pr=Ta/(94.72+170.75*G)
        return Ta,P_pr,T_pr

    #-----------------Pure Co2--------------------------#
    def rho_K_co2(P,T,G):
        """Compute CO2 properties as a function of temperature and pressure using modified Batzle-Wang equations

        Parameters
        ----------
        P : float or array-like
            Pressure in MPa
        T : float or array-like
            Temperature in °C
        G : float or array-like
            Gas gravity

        Returns
        -------
        float or array-like
            rho (g/cc): gas density 
            K (GPa): bulk modulus

        References
        ----------
            Xu, H. (2006). Calculation of CO2 acoustic properties using Batzle-Wang equations. Geophysics, 71(2), F21-F23.
        """        
           
        R=8.3145 # J.mol-1K-1 gas constant 
    
        Ta= T+273.15
        P_pr= P/7.4
        T_pr= Ta/(31.1+273.5)
    
        E=0.109*(3.85-T_pr)**2*np.exp(-(0.45+8*(0.56-1/T_pr)**2)*P_pr**1.2/T_pr)
        Z=(0.03+0.00527*(3.5-T_pr)**3)*P_pr+(0.642*T_pr-0.007*T_pr**4-0.52)+E 
        rho=28.8*G*P/(Z*R*Ta)

        r_0=0.85+5.6/(P_pr+2)+27.1/(P_pr+3.5)**2-8.7*np.exp(-0.65*(P_pr+1))
        #dz_dp=(0.03+0.00527*(3.5-T_pr)**3)+0.109*(3.85-T_pr)**2*1.2*P_pr**0.2*-(0.45+8*(0.56-1/T_pr)**2)/T_pr*np.exp(-(0.45+8*(0.56-1/T_pr)**2)*P_pr**1.2/T_pr)
        dzdp=BW.dz_dp(P_pr,T_pr)
        K=P/(1-P_pr*dzdp/Z)*r_0
        return rho, K/1000

    #-----------------Pure Gas--------------------------#
    def rho_K_gas(P, T, G):
        """Estimate the Gas density and bulk modulus at specific temperature and pressure. 

        Parameters
        ----------
        P : float or array-like
            Pressure in MPa
        T : float or array-like
            Temperature in °C
        G : float or array-like
            Gas gravity

        Returns
        -------
        float or array-like
            rho: Gas density
            K: Gas bulk modulus
        """        
   
        R=8.3145 # J.mol-1K-1 gas constant 
        Ta,P_pr,T_pr=BW.pseudo_p_t(P,T,G)
        E=0.109*(3.85-T_pr)**2*np.exp(-(0.45+8*(0.56-1/T_pr)**2)*P_pr**1.2/T_pr)
        Z=(0.03+0.00527*(3.5-T_pr)**3)*P_pr+(0.642*T_pr-0.007*T_pr**4-0.52)+E 
        rho=28.8*G*P/(Z*R*Ta)

        r_0=0.85+5.6/(P_pr+2)+27.1/(P_pr+3.5)**2-8.7*np.exp(-0.65*(P_pr+1))
        #dz_dp=(0.03+0.00527*(3.5-T_pr)**3)+0.109*(3.85-T_pr)**2*1.2*P_pr**0.2*-(0.45+8*(0.56-1/T_pr)**2)/T_pr*np.exp(-(0.45+8*(0.56-1/T_pr)**2)*P_pr**1.2/T_pr)
        dzdp=BW.dz_dp(P_pr,T_pr)
        K=P/(1-P_pr*dzdp/Z)*r_0
        return rho,K/1000

    #-----------------Pure Oil--------------------------#

    def rho_K_oil(P,T,den):
        """Estimate the oil density and bulk modulus at specific temperature and pressure.

        Parameters
        ----------
        P : float or array-like
            Pressure in MPa
        T : float or array-like
            Temperature in °C
        den : float or array-like
            oil density in g/cm3

        Returns
        -------
        float or array-like
            rho: oil density
            K: oil bulk modulus
        """        
            
        rho_p=den+(0.00277*P-1.71*0.0000001*P**3)*(den-1.15)**2+3.49*0.0001*P
        rho=rho_p/(0.972+3.81*0.0001*(T+17.78)**1.175)
        v=2096*(den/(2.6-den))**0.5-3.7*T+4.64*P+0.0115*(4.12*(1.08/den-1)**0.5-1)*T*P
        K=rho*v**2
        return rho, K

    #-----------------Gas satutrated with oil--------------------------#

    def rho_K_go(P, T, den,G,Rg):
        """compute density and bulk modulus of live oil. 

        Parameters
        ----------
        P : float or array-like
            Pressure in MPa
        T : float or array-like
            Temperature in °C
        den : float or array-like
            oil density in g/cm3
        G : float or array-like
            gas gravity 
        Rg : float or array-like
            the volume ratio of liberated gas to remaining oil at atmospheric pressure and 15.6°C, Liter/Liter

        Returns
        -------
        float or array-like
            v (m/s): velocity
            rho_g (g/cm3): true density of live oil at saturation
            K (GPa): true bulk modulus of live oil at saturation
        """          
        
        if Rg == None:
            Rg=0.02123*G*(P*np.exp(4.072/den-0.00377*T))**1.205
        
        B=0.972+0.00038*(2.4*Rg*(G/den)**0.5+T+17.8)**1.175
        rho_p=den*(1+0.001*Rg)**-1*B**-1 # pseudodensity
        v=2096*(rho_p/(2.6-rho_p))**0.5-3.7*T+4.64*P+0.0115*(4.12*(1.08/rho_p-1)**0.5-1)*T*P
        # true density of live oil at saturation
        rho_g=(den+0.0012*G*Rg)/B
        K=rho_g*v**2
        # return v
        return rho_g, K

    #-----------------Pure Water--------------------------#


    def rho_K_water(T, P):
        """Compute the density and bulk modulus of pure water as a function of temperature and pressure using Batzle and Wang (1992).

        Parameters
        ----------
        T : float or array-like
            Temperature in °C
        P : float or array-like
            Pressure in MPa

        Returns
        -------
        float or array-like
            rho_w (g/cm3): density of pure water
        """        
   
        #pressure=pressure*1e6
        # Align with the symbols and units in Batzle & Wang.
        #T, P = np.asanyarray(temperature), np.asanyarray(pressure)

        rho_w  = 1+1e-6*(-80*T- 3.3*T**2+0.00175*T**3+489*P-2*T*P+0.016*P*T**2-1.3e-5*T**3*P-0.333*P**2-0.002*T*P**2)
        v_w= BW.v_water(T, P)
        K_w= rho_w*v_w**2*1e-6
        return rho_w, K_w

    def v_water(T, P):
        """Acoustic velocity of pure water as a function of temperature
        and pressure using Batzle and Wang (1992).

        Parameters
        ----------
        T : float or array-like
            Temperature in °C
        P : float or array-like
            Pressure in MPa

        Returns
        -------
        float or array-like
            v_w (m/s): acoustic velocity of pure water 
        """        

        if np.any(P > 100):
            warnings.warn(
                'pressures above about 100 MPa-> inaccurate estimations'
            )
        w = np.array([[ 1.40285e+03,  1.52400e+00,  3.43700e-03, -1.19700e-05],
                    [ 4.87100e+00, -1.11000e-02,  1.73900e-04, -1.62800e-06],
                    [-4.78300e-02,  2.74700e-04, -2.13500e-06,  1.23700e-08],
                    [ 1.48700e-04, -6.50300e-07, -1.45500e-08,  1.32700e-10],
                    [-2.19700e-07,  7.98700e-10,  5.23000e-11, -4.61400e-13]])
        v_w=np.sum(w[i, j] * T**i * P**j for i in range(5) for j in range(4))
        
        return v_w

    def rho_K_brine(T,P,S):
        """Calculation of the density and bulk modulus of brine (NaCl) as a function of temperature, salinity and pressure using Batzle and Wang (1992).

        Parameters
        ----------
        T : float or array-like
            Temperature in °C
        P : float or array-like
            Pressure in MPa
        S : float or array-like
            weight fraction of sodium chloride in ppm/1e6

        Returns
        -------
        float or array-like
            rho_b (g/cm3): the density of brine
            K_b (GPa):bulk modulus of brine
        """        
   
        rho_w,_ = BW.rho_K_water(T, P)
        x = 300*P-2400*P*S+ T*(80+ 3*T- 3300*S - 13*P + 47*P*S)
        rho_b = rho_w + S*(0.668+ 0.44*S + 1e-6*x)
        v_b = BW.v_brine(T, P, S)
        K_b = rho_b*v_b**2*1e-6
        return rho_b, K_b

    def v_brine(T, P, S):
        """Calculte the acoustic velocity of brine as a function of temperature, salinity and pressure using Batzle and Wang (1992).

        Parameters
        ----------
        T : float or array-like
            Temperature in °C
        P : float or array-like
            Pressure in MPa
        S : float or array-like
            weight fraction of sodium chloride in ppm/1e6

        Returns
        -------
        float or array-like
            v_b (m/s): the velocity of brine
        """        
         
        v_w = BW.v_water(T, P)
        s1 = 1170 - 9.6*T + 0.055*T**2 - 8.5e-5*T**3 + 2.6*P - 0.0029*T*P - 0.0476*P**2
        s15 = 780 - 10*P + 0.16*P**2
        s2 = -820
        v_b = v_w + s1 * S + s15 * S**1.5 + s2 * S**2
        return v_b


