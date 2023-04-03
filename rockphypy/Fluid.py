import scipy.special as sp
import numpy as np
from rockphypy.utils import utils
#from utils import *


class Fluid:
    """
    Fluid subsitution approaches and models describing velocity dispersion and attenuation due to the fluid effect. 
    """ 

    @staticmethod 
    def Brie(Kw, Kgas, Sw, e):
        """Brie empirical fluid mixing law

        Parameters
        ----------
        Kw : float
            bulk modulus of fluid phase
        Kgas : float
            bulk modulus of gas phase
        Sw : float or array
            water saturation
        e : int
            Brie component

        Returns
        -------
        float or array
            Kf: effective fluid propertie
        """        
        
        Kf= (Kw-Kgas)*Sw**e+Kgas 
        return Kf   

    @staticmethod    
    def Biot(Kdry,Gdry,K0,Kfl,rho0,rhofl,eta,phi,kapa,a,alpha,freq):
        """
        Compute Biot dispersion and velocity attenuation 

        Parameters
        ----------
        Kdry : float or array-like
            dry frame bulk modulus
        Gdry : float or array-like
            dry frame shear modulus
        K0 : float
            bulk modulus of mineral material making up rock
        Kfl : float
            effective bulk modulus of pore fluid
        rho0 : float
            grain density
        rhofl : float
            pore fluid density
        eta : float
            η is the viscosity of the pore fluid
        phi : float
            porosity 
        kapa : float
            absolute permeability of the rock
        a : float
            pore-size parameter. Stoll (1974) found that values between 1/6 and 1/7 of the mean grain diameter
        alpha : float
            tortuosity parameter, always greater than or equal to 1.
        freq : float or array-like
            frequency range, e.g 10^-3 to 10^3 Hz

        Returns
        -------
        float or array-like
            Vp_fast, fast P-wave velocities at all frequencies
            Vp_slow, slow P-wave velocities at all frequencies
            Vs, S-wave velocities
            QP1_inv, fast P-wave attenuation
            QP2_inv, slow P-wave attenuation
            Qs_inv, S-wave attenuation
        """        

        rho=(1-phi)*rho0+phi*rhofl # bulk density
        # poroelastic coefficients
        D = K0*(1+phi*(K0/Kfl-1))
        M = K0**2/(D-Kdry)
        C = (K0-Kdry)*K0/(D-Kdry)
        H = Kdry+4*Gdry/3+(K0-Kdry)**2/(D-Kdry)
        # angular frequency of the plane wave, ω is the angular frequency of the plane wave.
        w = 2*np.pi*freq 
        # viscodynamic operator
        zeta = np.sqrt(w*a**2*rhofl/eta)
        T = np.zeros(zeta.size, dtype=complex)
        # asymptotic values for very high and very low frequencies
        # high frequency limit, T->(1+i)/sqrt(2)
        T[zeta<=1e3]=np.exp(1j*3*np.pi/4)*sp.jv(1,zeta[zeta<=1e3]*np.exp(-1j*np.pi/4) )/sp.jv(0, zeta[zeta<=1e3]*np.exp(-1j*np.pi/4) )
        T[zeta>1e3]=(1+1j)/np.sqrt(2)
        # low frequency limit F->1
        F = 1/4 * (zeta*T/(1+2*1j*T/zeta))
        F[zeta<1e-1]=1
        q = alpha*rhofl/phi - 1j*eta*F/(w*kapa)
        # a,b,c terms of the root 
        Ta = (C**2-M*H)
        Tb = H*q+M*rho-2*C*rhofl
        Tc = rhofl**2-rho*q
        # slowness squared
        P1_slowness2= (-Tb+np.sqrt(Tb**2-4*Ta*Tc) )/(2*Ta)
        P2_slowness2= (-Tb-np.sqrt(Tb**2-4*Ta*Tc) )/(2*Ta)
        S_slowness2= (rho*q-rhofl**2)/(Gdry*q)
        # velocities and attenuations
        Vp_fast = 1/ np.sqrt(P1_slowness2).real
        Vp_slow = 1/ np.sqrt(P2_slowness2).real
        Vs = 1/ np.sqrt(S_slowness2).real
        QP1_inv= (1/P1_slowness2).imag/(1/P1_slowness2).real
        QP2_inv= (1/P2_slowness2).imag/(1/P2_slowness2).real
        Qs_inv = (1/S_slowness2).imag/(1/S_slowness2).real

        return Vp_fast,Vp_slow,Vs,QP1_inv,QP2_inv,Qs_inv

    @staticmethod        
    def Biot_HF(Kdry,Gdry,K0,Kfl,rho0,rhofl,phi,alpha):
        """Biot high-frequency limiting velocities in the notation of Johnson and Plona (1982)

        Parameters
        ----------
        Kdry : float or array-like
            dry frame bulk modulus
        Gdry : float or array-like
            dry frame shear modulus
        K0 : float
            bulk modulus of mineral material making up rock 
        Kfl : float
            effective bulk modulus of pore fluid
        rho0 : float
            grain density
        rhofl : float
            pore fluid density
        phi : float
            porosity 
        alpha : float
            tortuosity parameter, always greater than or equal to 1.

        Returns
        -------
        float or array-like
            Vp_fast,Vp_slow,Vs:  high-frequency limiting velocities
        """        
      
        rho=(1-phi)*rho0+phi*rhofl # bulk density
        rho12=(1-alpha)*phi*rhofl
        rho22=alpha*phi*rhofl
        rho11=(1-phi)*rho0-(1-alpha)*phi*rhofl
        # repetative terms
        T1=1-phi-Kdry/K0
        T2= phi*K0/Kfl
        R=phi**2*K0/(T1+T2)
        Q=(T1*phi*K0/(T1+T2))
        P=((1-phi)*T1*K0+T2*Kdry)/(T1+T2) + 4*Gdry/3
        Delta= P*rho22+R*rho11-2*Q*rho12
        T3= rho11*rho22-rho12**2
        T4=P*R-Q**2
        Vp_fast= np.sqrt( (Delta+np.sqrt(Delta**2-4*T3*T4))/(2*T3))
        Vp_slow= np.sqrt( (Delta-np.sqrt(Delta**2-4*T3*T4))/(2*T3))
        Vs=np.sqrt(Gdry/(rho-phi*rhofl*alpha**-1))
        return Vp_fast,Vp_slow,Vs

    @staticmethod
    def Geertsma_Smit_HF(Kdry,Gdry,K0,Kfl,rho0,rhofl,phi,alpha):
        """Approximation of Biot high-frequency limit of the fast P-wave velocity given by Geertsma and Smit (1961), This form predicts velocities that are too high (by about 3%–6%) compared with the actual high-frequency limit.

        Parameters
        ----------
        Kdry : float or array-like
            dry frame bulk modulus
        Gdry : float or array-like
            dry frame shear modulus
        K0 : float
            bulk modulus of mineral material making up rock 
        Kfl : float
            effective bulk modulus of pore fluid
        rho0 : float
            grain density
        rhofl : float
            pore fluid density
        phi : float
            porosity 
        alpha : float
            tortuosity parameter, always greater than or equal to 1.

        Returns
        -------
        float or array-like
            Vp_fast,Vs: high-frequency limiting velocities
        """        
           
        rho=(1-phi)*rho0+phi*rhofl # bulk density
        rho_biot= rho0*(1-phi)+phi*rhofl*(1-alpha**-1)
        Hdry= Kdry+4*Gdry/3
        T1=phi*rho/(rhofl*alpha)
        alpha_biot= 1-Kdry/K0
        Vp_fast=np.sqrt( 1/rho_biot*(Hdry+(T1+alpha_biot*(alpha_biot-2*phi/alpha))/((alpha_biot-phi)/K0+phi/Kfl )))

        Vs=np.sqrt(Gdry/rho_biot)

        return Vp_fast,Vs

    @staticmethod
    def Geertsma_Smit_LF(Vp0,Vpinf, freq,phi, rhofl, kapa, eta):
        """Low and middle-frequency approximations of Biot wave given by Geertsma and Smit (1961). Noticed that mathematically this approximation is valid at moderate-to-low seismic frequencies, i.e. f<fc

        Parameters
        ----------
        Vp0 : float
            Biot−Gassmann low-frequency limiting P-wave velocity
        Vpinf : float
            Biot highfrequency limiting P-wave velocity
        freq : float or array-like
            frequency
        phi : float
            porosity
        rhofl : float
            fluid density
        kapa : float
            absolute permeability of the rock.
        eta : float
            viscosity of the pore fluid

        Returns
        -------
        float or array-like
            Vp: frequency-dependent P-wave velocity of saturated rock
        """        
           
        # Biot reference frequency
        fc = phi*eta/(2*np.pi*rhofl*kapa)
        # if freq >=fc:
        #     print('ATTENTION, the input frequency is higher than the reference frequency, the approximation might not be valid!')
        a = (fc/freq)**2
        Vp= np.sqrt((Vpinf**4+Vp0**4*a)/(Vpinf**2+Vp0**2*a))
        return Vp

    @staticmethod
    def Gassmann(K_dry,G_dry,K_mat,Kf,phi):
        """Computes saturated elastic moduli of rock via Gassmann equation given dry-rock moduli. 

        Parameters
        ----------
        K_dry : float or array-like
            dry frame bulk modulus 
        G_dry : float or array-like
            dry frame shear modulus 
        K_mat : float
            matrix bulk modulus
        Kf : float
            fluid bulk modulus
        phi : float or array-like
            porosity

        Returns
        -------
        float or array-like
            K_sat, G_sat: fluid saturated elastic moduli
        """        
          
        A=(1-K_dry/K_mat)**2
        B=phi/Kf+(1-phi)/K_mat-K_dry/(K_mat**2)
        K_sat=K_dry+A/B
        G_sat = G_dry # At low frequencies, Gassmann’s relations predict no change in the shear modulus between dry and saturated patches
        return K_sat,G_sat

    @staticmethod
    def Gassmann_sub(phi, K0, Ksat1,Kfl1,Kfl2):
        """Fluid subsititution using Gassmann equation, thr rock is initially saturated with a fluid, compute the saturated moduli for tge rock saturated with a different fluid

        Parameters
        ----------
        phi : float or array-like 
            porosity
        K0 : float 
            mineral modulus
        Ksat1 : float or array-like
            original bulk modulus of rock saturated with fluid of bulk modulus Kfl1
        Kfl1 : float 
            original saturant
        Kfl2 : float 
            new saturant

        Returns
        -------
        float or array-like
            Ksat2: new satuarted bulk modulus of the rock 
        """       

        a=Ksat1/(K0-Ksat1)-Kfl1/(phi*(K0-Kfl1))+Kfl2/(phi*(K0-Kfl2))
        Ksat2= a*K0/(1+a)

        return Ksat2

    @staticmethod
    def vels(K_dry,G_dry,K0,D0,Kf,Df,phi):
        """Computes Vp,Vs and densities of saturated rock using Gassmann relations from elastic moduli of rock. See also `Gassmann_vels`.

        Parameters
        ----------
        K_dry : float
            dry frame bulk modulus
        G_dry : float
            dry frame shear modulus
        K0 : float
            mineral matrix bulk modulus
        D0 : float
            mineral matrix density
        Kf : float
            fluid bulk modulus
        Df : float
            fluid density in g/cm3
        phi : float or array
            porosity

        Returns
        -------
        float or array
            Vp, Vs, rho
        """        
         
        #D0=Dsh*vsh+Dqz*(1-vsh)####
        rho  = D0*(1-phi)+Df*phi
        K,_= Fluid.Gassmann(K_dry,G_dry,K0,Kf,phi)
        Vp   = np.sqrt((K+4./3*G_dry)/rho)*1e3
        Vs   = np.sqrt(G_dry/rho)*1e3
        return Vp, Vs, rho
    @staticmethod
    def Gassmann_vels(Vp1,Vs1,rho1,rhofl1,Kfl1,rhofl2,Kfl2,K0,phi):
        """Gassmann fluid substituion with velocities 

        Parameters
        ----------
        Vp1 : float or array-like
            saturated P velocity of rock with fluid 1
        Vs1 : float or array-like
            saturated S velocity of rock with fluid 1
        rho1 : float
            bulk density of saturated rock with fluid 1
        rhofl1 : float
            density of fluid 1
        Kfl1 : float
            bulk modulus of fluid 1
        rhofl2 : float
            density of fluid 2
        Kfl2 : float
            bulk modulus of fluid 2
        K0 : float
            mineral bulk modulus
        phi : float or array-like
            porosity

        Returns
        -------
        float or array-like
            Vp2, Vs2: velocities of rock saturated with fluid 2
        """        
            
        rho2=rho1 - phi*rhofl1 +phi*rhofl2 # density update
        G1=rho1*Vs1**2 
        Ksat1=rho1*Vp1**2-(4/3)*G1
        k2= Fluid.Gassmann_sub(phi, K0, Ksat1,Kfl1,Kfl2)
        G2=G1
        Vp2=np.sqrt((k2+(4/3)*G2)/rho2) 
        Vs2=np.sqrt(G2/rho2)
        return Vp2, Vs2

    @staticmethod
    def Gassmann_approx(Msat1,M0,Mfl1,phi,Mfl2):
        """Perform gassmann fluid subsititution using on p wave modulus 

        Parameters
        ----------
        Msat1 : float or array-like
            in situ saturated p wave modulus from well log data
        M0 : float
            mineral modulus
        Mfl1 : float
            p wave modulus of in situ fluid 
        phi : float
            porosity
        Mfl2 : float
            p wave modulus of new fluid for susbtitution

        Returns
        -------
        float or array-like
            Msat2:  p wave modulus of rock fully saturated with new fluid
        """        
            
        x=Msat1/(M0-Msat1) -Mfl1/(phi*(M0-Mfl1)) +Mfl2/(phi*(M0-Mfl2))
        Msat2=x*M0/(1+x)
        return Msat2

    @staticmethod
    def Brown_Korringa_dry2sat(Sdry,K0,G0,Kfl,phi):
        """Compute fluid saturated compliances from dry compliance for anisotropic rock using Brown and Korringa (1975). See eq. 32 in the paper.

        Parameters
        ----------
        Sdry : 2d array
            comliance matrix of the dry rock
        K0 : float
            Isotropic mineral bulk modulus
        G0 : float
            Isotropic mineral shear modulus
        Kfl : float
            Isotropic fluid bulk modulus
        phi : float
            porosity

        Returns
        -------
        2d array
            Ssat (6x6 matrix): Saturated compliance of anisotropic rock
        """        
           
        # Compressibilities of the fluid, mineral, and dry rock, sum over repeted index, e,g, β0 =s0[ααγγ]
        beta0 = 1/K0
        
        betadry = np.sum(Sdry[0:3,0:3]) 
        betafl =1/Kfl
        # Isotropic mineral compliance tensor
        C = utils.write_iso(K0,G0)
        S0 = np.linalg.inv(C)
        # nominator
        Sprime = np.sum(Sdry[0:3,:],axis=0) - np.sum(S0[0:3,:],axis=0)
        Sprime = Sprime.reshape((1,6))
        # denominator
        denom = (betafl-beta0)*phi+(betadry-beta0)
        Ssat=Sdry-np.dot(Sprime.T,Sprime)/denom
        return Ssat

    @staticmethod
    def Brown_Korringa_sat2dry(Ssat,K0,G0,Kfl,phi):
        """Compute dry compliance from fluid saturated compliances for arbitrarily anisotropic rock using Brown and Korringa (1975). See eq. 32 in the paper. 

        Parameters
        ----------
        Ssat : 2d array
            comliance matrix (6x6) of the saturated rock 
        K0 : float
            Isotropic mineral bulk modulus
        G0 : float
            Isotropic mineral shear modulus
        Kfl : float
            Isotropic fluid bulk modulus
        phi : float
            porosity

        Returns
        -------
        2d array
            Sdry (6x6 matrix): Dry compliance of anisotropic rock
        """        
      
        # Compressibilities of the fluid, mineral, and dry rock, sum over repeted index, e,g, β0 =s0[ααγγ]
        beta0 = 1/K0
        
        betasat = np.sum(Ssat[0:3,0:3]) 
        betafl =1/Kfl
        # Isotropic mineral compliance tensor
        C = utils.write_iso(K0,G0)
        S0 = np.linalg.inv(C)
        # nominator
        Sprime = np.sum(Ssat[0:3,:],axis=0) - np.sum(S0[0:3,:],axis=0)
        Sprime = Sprime.reshape((1,6))
        # denominator
        denom = (betafl-beta0)*phi-(betasat-beta0)
        Sdry=Ssat + np.dot(Sprime.T,Sprime)/denom
        return Sdry

    @staticmethod
    def Brown_Korringa_sub(Csat,K0,G0,Kfl1,Kfl2,phi):
        """Fluid substitution in arbitrarily anisotropic rock using Brown and Korringa (1975). the rock is originally saturated by fluid 1. After fluid subsititution, the rock is finally saturated by fluid 2.

        Parameters
        ----------
        Csat : 6x6 matrix
            comliance matrix of the saturated rock
        K0 : float
            Isotropic mineral bulk modulus
        G0 : float
            Isotropic mineral shear modulus
        Kfl1 : float
            bulk modulus of the original fluid
        Kfl2 : float
            bulk modulus of the final fluid
        phi : float
            porosity

        Returns
        -------
        2d array
            Csat2, Ssat2 (6x6 matrix): Dry stiffness and compliance matrix of anisotropic rock saturated with new fluid
        """        
        Ssat = np.linalg.inv(Csat)
        
        Sdry = Fluid.Brown_Korringa_sat2dry(Ssat,K0,G0,Kfl1,phi)

        Ssat2 = Fluid.Brown_Korringa_dry2sat(Sdry,K0,G0,Kfl2,phi)
        Csat2 = np.linalg.inv(Ssat2)
        
        return Csat2, Ssat2

    @staticmethod
    def Mavko_Jizba(Vp_hs, Vs_hs,Vpdry, Vsdry, K0, rhodry, rhofl,Kfl, phi):
        """Predicting the very high-frequency moduli and velocities of saturated rocks from dry rock properties using the Squirt flow model derived by Mavko and Jizba (1991). 

        Parameters
        ----------
        Vp_hs : float
            P wave velocity of the dry rock measured at very high effective pressure in the unit of m/s
        Vs_hs : float
            S wave velocity of the dry rock  measured at very high effective pressure in the unit of m/s
        Vpdry : array
            P wave velocity of the dry rock measured at different effective pressure in the unit of m/s
        Vsdry : array
            S wave velocity of the dry rock measured at different effective pressure in the unit of m/s
        K0 : float
            mineral bulk moduli
        rhodry : float
            bulk density of the dry rock
        rhofl : float
            bulk density of the pore fluid
        Kfl : float
            bulk moduli of the pore fluid
        phi : float
            porosity

        Returns
        -------
        _type_
            Kuf_sat (float):GPa, predicted high frequency bulk moduli of saturated rock
            Guf_sat (array): GPa, predicted high frequency shear moduli of saturated rock at different pressure 
            Vp_hf (array): m/s, predicted high frequency P wave velocities of saturated rock
            Vs_hf (array): m/s, predicted high frequency S wave velocities of saturated rock
        """        
        # dry rock moduli with pressure 
        Kdry,Gdry = utils.M_from_V(rhodry, Vpdry,Vsdry)
        # dry moduli dry rock at high high effective pressure, crack free high pressure moduli
        Khs, Ghs = utils.M_from_V(rhodry, Vp_hs,Vs_hs)
        # high-frequency “wet-frame moduli
        Kuf= Khs # ignore the a first order correction for the difference in fluid and mineral compressibilities. 
        # high-frequency saturated moduli 
        Kuf_sat,_ =  Fluid.Gassmann(Kuf,1,K0,Kfl,phi)
        Guf_sat_inv = 1/Gdry - 4/15*(1/Kdry-1/Kuf)
        Guf_sat = 1/Guf_sat_inv
        # predicted high frequency velocities
        rho_sat=rhodry+phi*rhofl # saturated bulk density
        Vp_hf,Vs_hf = utils.V(Kuf_sat, Guf_sat, rho_sat)
        return Kuf_sat,Guf_sat,Vp_hf,Vs_hf

    @staticmethod
    def Squirt_anisotropic(Sdry, Sdry_hp):
        """Predict wet unrelaxed frame compliances at very high frequency from dry frame compliances for transversely isotropic rocks using theoretical formula derived by Mukerji and Mavko, (1994)

        Parameters
        ----------
        Sdry : array
            dry rock compliances [S11 S12 S13 S33 S44]
        Sdry_hp : array
            dry rock compliances at very high effective stress [S11 S12 S13 S33 S44]

        Returns
        -------
        array
            The wet-frame compliances [S11 S12 S13 S33 S44]
        """        
        Sdry=np.asanyarray(Sdry)
        Sdry_hp=np.asanyarray(Sdry_hp)
        # delta_Sdry_ijkl
        DSdry=Sdry-Sdry_hp
        # delta_Sdry_aabb
        DSdry_aabb= 2*DSdry[0]+2*DSdry[1]+4*DSdry[2]+DSdry[3]
        #2*(dss(:,1)+dss(:,2)+2*dss(:,3))+dss(:,4);
        # Delta S_tilde 
        D_Sdry= DSdry/DSdry_aabb
        # delta_Sdry_abab= 2S11+S33+4S44+2S66
        D_Sdry_abab=2*D_Sdry[0]+D_Sdry[3]+4*D_Sdry[4]+4*(D_Sdry[0]-D_Sdry[1])

        alpha=0.25*(D_Sdry_abab-1)

        #the elements of Gijkl
        T = 1-4*alpha
        G11=D_Sdry[0]-4*alpha/T*(D_Sdry[1]+D_Sdry[2])
        G12=D_Sdry[1]/T 
        G13=D_Sdry[2]/T
        G33=D_Sdry[3]-8*alpha*D_Sdry[2]/T
        G44=D_Sdry[4]/T-(D_Sdry[0]+D_Sdry[2])/(4*T)+(G11+G33)/4
        # G66=2*(G11-G12)

        G=np.array([G11,G12,G13,G33,G44])

        # The wet-frame compliance 
        S_wet= Sdry-DSdry_aabb*G
        return S_wet

    @staticmethod
    def White_Dutta_Ode(Kdry, Gdry, K0, phi, rho0, rhofl1,rhofl2, Kfl1, Kfl2,eta1,eta2,kapa,a,sg,freq ):
        """Dispersion and Attenuation of partial saturation using White and Dutta–Odé Model. 

        Parameters
        ----------
        Kdry : float
            bulk modulus of the dry rock 
        Gdry : float
            shear modulus of the dry rock 
        K0 : float
            Isotropic mineral bulk modulus
        phi : float
            porosity
        rho0 : float
            mineral density
        rhofl1 : float
            density of the fluid opcupying the central sphere
        rhofl2 : float
            density of the fluid opcupying the outer sphere
        Kfl1 : float
            bulk modulus of the fluid opcupying the central sphere
        Kfl2 : float
            bulk modulus of the fluid opcupying the outer sphere
        eta1 : float
            viscousity of the fluid opcupying the central sphere
        eta2 : float
            viscousity of the fluid opcupying the outer sphere
        kapa : float
            absolute permeability of the rock
        a : float
            radius of central sphere , sg=a3/b3
        sg : float
            saturation of fluid opcupying the central sphere
        freq : float or array-like
            frequencies 

        Returns
        -------
        float, array-like
            Vp (m/s): P wave velocity 
            a_w: attenuation coefficient
            K_star: complex bulk modulus 
        """        
        
        omega= 2*np.pi*freq
        #sg=a3/b3, outer radius b
        b= a/sg**(1/3)
        # The saturated bulk and shear moduli, Kj and μj of region j
        K1,_=Fluid.Gassmann(Kdry,0,K0,Kfl1,phi)
        K2,_=Fluid.Gassmann(Kdry,0,K0,Kfl2,phi)
        G1=Gdry
        G2=G1
        # 
        R1=(K1-Kdry)/(1-Kdry/K0) * (3*K2+4*G2)/(K2*(3*K1+4*G2)+4*G2*(K1-K2)*sg )
        R2=(K2-Kdry)/(1-Kdry/K0) * (3*K1+4*G1)/(K2*(3*K1+4*G2)+4*G2*(K1-K2)*sg )
        KA1=(phi/Kfl1+(1-phi)/K0-Kdry/K0**2)**-1
        KA2=(phi/Kfl2+(1-phi)/K0-Kdry/K0**2)**-1
        KE1=(1-Kfl1*(1-K1/K0)*(1-Kdry/K0)/(phi*K1*(1-Kfl1/K0)))*KA1
        KE2=(1-Kfl2*(1-K2/K0)*(1-Kdry/K0)/(phi*K2*(1-Kfl2/K0)))*KA2
        alpha1=(1j*omega*eta1/(kapa*KE1))**(1/2)
        alpha2=(1j*omega*eta2/(kapa*KE2))**(1/2)

        Z1=eta1*a/kapa *(1-np.exp(-2*alpha1*a))/((alpha1*a-1)+(alpha1*a+1)*np.exp(-2*alpha1*a)   )
        Z2=-eta2*a/kapa *(   (alpha2*b+1) +(alpha2*b-1)*np.exp(2*alpha2*(b-a) ) )/(  (alpha2*b+1)* (alpha2*a-1)- (alpha2*b-1)*(alpha2*a+1)*np.exp(2*alpha2*(b-a))  )
        Q1= (1-Kdry/K0)*KA1/K1
        Q2= (1-Kdry/K0)*KA2/K2
        W=3*a**2*(R1-R2)*(-Q1+Q2)/ (b**3*1j*omega*(Z1+Z2) )
        Kinf= (K2*(3*K1+4*G2)+4*G2*(K1-K2)*sg)/((3*K1+4*G2)-3*(K1-K2)*sg)
        Klowf = (K2*(K1-Kdry)+sg*Kdry*(K2-K1))/((K1-Kdry)+sg*(K2-K1))

        K_star = Kinf/(1-Kinf*W)

        # when the central sphere is saturated with a very compressible gas
        K_star_gas= Kinf/(1+3*a**2*R2*Q2*Kinf/(b**3*1j*omega*Z2) )

        # P wave modulus and density 
        M= K_star+4*Gdry/3
        rho = (1-phi)*rho0+phi*sg*rhofl1+phi*(1-sg)*rhofl2
        # phase angle 
        theta= np.arctan(M.imag/M.real)
        # velocity and attenuation coefficient, see 3.8.56 in RPH 
        Vp = np.sqrt(np.abs(M)/rho) * 1/np.cos(theta/2)
        a_w= omega/Vp *np.tan(theta/2)

        return Vp,a_w, K_star
