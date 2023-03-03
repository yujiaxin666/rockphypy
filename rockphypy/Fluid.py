import scipy.special as sp
import numpy as np
from rockphypy.utils import utils
#from utils import *
class Fluid(utils):
    """
    
    """    
    def Biot(Kdry,Gdry,K0,Kfl,rho0,rhofl,eta,phi,kapa,a,alpha,freq):
        """_summary_

        Args:
            Kdry (scalar or array): dry frame bulk modulus
            Gdry (scalar or array): dry frame shear modulus
            K0 (scalar or array):   bulk modulus of mineral material making up rock
            Kfl (GPa): effective bulk modulus of pore fluid
            rho0 (g/cm3): grain density
            rhofl (g/cm3): pore fluid density
            eta : η is the viscosity of the pore fluid
            phi (frac): porosity 
            kapa : κ is the absolute permeability of the rock
            a : pore-size parameter. Stoll (1974) found that values between 1/6 and 1/7 of the mean grain diameter
            alpha: tortuosity parameter, always greater than or equal to 1.

            freq (scalar or array): frequency range, e.g 10^-3 to 10^3 Hz


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
        
    def Biot_HF(Kdry,Gdry,K0,Kfl,rho0,rhofl,phi,alpha):
        """Biot high-frequency limiting velocities in the notation of Johnson and Plona (1982)

        Args:
            Kdry (scalar or array): dry frame bulk modulus
            Gdry (scalar or array): dry frame shear modulus
            K0 (scalar or array):   bulk modulus of mineral material making up rock
            Kfl (GPa): effective bulk modulus of pore fluid
            rho0 (g/cm3): grain density
            rhofl (g/cm3): pore fluid density
            phi (frac): porosity 
            alpha: tortuosity parameter, always greater than or equal to 1.

        Returns:
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

    def Geertsma_Smit_HF(Kdry,Gdry,K0,Kfl,rho0,rhofl,phi,alpha):
        """ Approximation of Biot high-frequency limit of the fast P-wave velocity given by Geertsma and Smit (1961), This form predicts velocities that are too high (by about 3%–6%) compared with the actual high-frequency limit.

        Args:
            Kdry (scalar or array): dry frame bulk modulus
            Gdry (scalar or array): dry frame shear modulus
            K0 (scalar or array):   bulk modulus of mineral material making up rock
            Kfl (GPa): effective bulk modulus of pore fluid
            rho0 (g/cm3): grain density
            rhofl (g/cm3): pore fluid density
            phi (frac): porosity 
            alpha: tortuosity parameter, always greater than or equal to 1.

        Returns:
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

    def Geertsma_Smit_LF(Vp0,Vpinf, freq,phi, rhofl, kapa, eta):
        """Low and middle-frequency approximations of Biot wave given by Geertsma and Smit (1961). Noticed that mathematically this approximation is valid at moderate-to-low seismic frequencies, i.e. f<fc

        Args:
            Vp0 (_type_): Biot−Gassmann low-frequency limiting P-wave velocity
            Vpinf (_type_): Biot highfrequency limiting P-wave velocity
            freq (_type_): frequency
            phi (_type_): porosity
            rhofl (_type_): fluid density
            kapa (_type_): absolute permeability of the rock.
            eta (_type_): viscosity of the pore fluid

        Returns:
            Vp: frequency-dependent P-wave velocity of saturated rock
        """    
        # Biot reference frequency
        fc = phi*eta/(2*np.pi*rhofl*kapa)
        # if freq >=fc:
        #     print('ATTENTION, the input frequency is higher than the reference frequency, the approximation might not be valid!')
        a = (fc/freq)**2
        Vp= np.sqrt((Vpinf**4+Vp0**4*a)/(Vpinf**2+Vp0**2*a))
        return Vp

    def Gassmann(K_dry,G_dry,K_mat,Kf,phi):
        """Computes saturated elastic moduli of rock via Gassmann equation given dry-rock moduli. 

        Args:
            K_dry (Gpa): dry frame bulk modulus 
            G_dry (Gpa): dry frame shear modulus 
            K_mat (Gpa): matrix bulk modulus
            Kf (Gpa): fluid bulk modulus
            phi (frac): porosity

        Returns:
            K_sat, G_sat: fluid saturated elastic moduli
        Written by Jiaxin Yu
        """    
        A=(1-K_dry/K_mat)**2
        B=phi/Kf+(1-phi)/K_mat-K_dry/(K_mat**2)
        K_sat=K_dry+A/B
        G_sat = G_dry # At low frequencies, Gassmann’s relations predict no change in the shear modulus between dry and saturated patches
        return K_sat,G_sat

    def Gassmann_sub(phi, K0, Ksat1,Kfl1,Kfl2):
        """Fluid subsititution using Gassmann equation, thr rock is initially saturated with a fluid, compute the saturated moduli for tge rock saturated with a different fluid
        Args:
            phi (frac): porosity
            K0 (GPa): mineral modulus
            Ksat1 (GPa): original bulk modulus of rock saturated with fluid of bulk modulus Kfl1
            Kfl1 (GPa): original saturant
            Kfl2 (GPa): new saturant

        Returns:
            Ksat2: new satuarted bulk modulus of the rock 
        """    

        a=Ksat1/(K0-Ksat1)-Kfl1/(phi*(K0-Kfl1))+Kfl2/(phi*(K0-Kfl2))
        Ksat2= a*K0/(1+a)

        return Ksat2
    def Gassmann_vels(Vp1,Vs1,rho1,rhofl1,Kfl1,rhofl2,Kfl2,K0,phi):
        """Gassmann fluid substituion with velocities 

        Args:
            Vp1 (_type_): saturated P velocity of rock with fluid 1
            Vs1 (_type_): saturated S velocity of rock with fluid 1
            rho1 (_type_): bulk density of saturated rock with fluid 1
            rhofl1 (_type_):density of fluid 1
            Kfl1 (_type_): bulk modulus of fluid 1
            rhofl2 (_type_):density of fluid 2
            Kfl2 (_type_): bulk modulus of fluid 2
            K0 (_type_): mineral bulk modulus
            phi (_type_): porosity

        Returns:
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


    def Gassmann_approx(Msat1,M0,Mfl1,phi,Mfl2):
        """ Perform gassmann fluid subsititution using on p wave modulus 

        Args:
            Msat1 (GPa): in situ p wave modulus from well log
            M0 (GPa): mineral modulus
            Mfl1 (GPa): p wave modulus of in situ fluid 
            phi (frac): por
            Mfl2 (GPa): p wave modulus of new fluid after susbtitution

        Returns:
            Msat2:  p wave modulus of fully saturated rock
        """    
        x=Msat1/(M0-Msat1) -Mfl1/(phi*(M0-Mfl1)) +Mfl2/(phi*(M0-Mfl2))
        Msat2=x*M0/(1+x)
        return Msat2

    def Brown_Korringa_dry2sat(Sdry,K0,G0,Kfl,phi):
        """Compute fluid saturated compliances from dry compliance for anisotropic rock using Brown and Korringa (1975). See eq. 32 in the paper. 

        Args:
            Sdry (6x6 matrix): comliance matrix of the dry rock
            K0 (GPa): Isotropic mineral bulk modulus
            G0 (GPa): Isotropic mineral shear modulus
            Kfl (GPa): Isotropic fluid bulk modulus
            phi (frac): porosity
        Returns:
            Ssat (6x6 matrix): Saturated compliance of anisotropic rock
        Written by Jiaxin Yu
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

    def Brown_Korringa_sat2dry(Ssat,K0,G0,Kfl,phi):
        """ Compute dry compliance from fluid saturated compliances for arbitrarily anisotropic rock using Brown and Korringa (1975). See eq. 32 in the paper. 

        Args:
            Ssat (6x6 matrix): comliance matrix of the saturated rock
            K0 (GPa): Isotropic mineral bulk modulus
            G0 (GPa): Isotropic mineral shear modulus
            Kfl (GPa): Isotropic fluid bulk modulus
            phi (frac): porosity
        Returns:
            Sdry (6x6 matrix): Dry compliance of anisotropic rock
        Written by Jiaxin Yu
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

    def Brown_Korringa_sub(Csat,K0,G0,Kfl1,Kfl2,phi):
        """ Fluid substitution in arbitrarily anisotropic rock using Brown and Korringa (1975). 
        the rock is originally saturated by fluid 1. After fluid subsititution, the rock is finally saturated by fluid 2.
        Args:
            Ssat (6x6 matrix): comliance matrix of the saturated rock
            K0 (GPa): Isotropic mineral bulk modulus
            G0 (GPa): Isotropic mineral shear modulus
            Kfl1 (GPa): bulk modulus of the original fluid
            Kfl2 (GPa): bulk modulus of the final fluid
            phi (frac): porosity
        Returns:
            Csat2, Ssat2 (6x6 matrix): Dry stiffness and compliance matrix of anisotropic rock saturated with new fluid
        Written by Jiaxin Yu
        """    
        Ssat = np.linalg.inv(Csat)
        
        Sdry = Fluid.Brown_Korringa_sat2dry(Ssat,K0,G0,Kfl1,phi)

        Ssat2 = Fluid.Brown_Korringa_dry2sat(Sdry,K0,G0,Kfl2,phi)
        Csat2 = np.linalg.inv(Ssat2)
        
        return Csat2, Ssat2

    def Mavko_Jizba(Vp_hs, Vs_hs,Vpdry, Vsdry, K0, rhodry, rhofl,Kfl, phi):
        """Predicting the very high-frequency moduli and velocities of saturated rocks from dry rock properties using the Squirt flow model derived by Mavko and Jizba (1991). 

        Args:
            Vp_hs (scalar): P wave velocity of the dry rock measured at very high effective pressure in the unit of m/s
            Vs_hs (scalar): S wave velocity of the dry rock  measured at very high effective pressure in the unit of m/s
            Vpdry (array): P wave velocity of the dry rock measured at different effective pressure in the unit of m/s
            Vsdry (array): S wave velocity of the dry rock measured at different effective pressure in the unit of m/s
            K0 (GPa): mineral bulk moduli
            rhodry (g/cm3): bulk density of the dry rock
            rhofl (g/cm3): bulk density of the pore fluid
            Kfl (GPa): bulk moduli of the pore fluid
            phi (frac): porosity

        Returns:
            Kuf_sat (scalar):GPa, predicted high frequency bulk moduli of saturated rock
            Guf_sat (array): GPa, predicted high frequency shear moduli of saturated rock at different pressure 
            Vp_hf (array): m/s, predicted high frequency P wave velocities of saturated rock
            Vs_hf (array): m/s, predicted high frequency S wave velocities of saturated rock
        Written by Jiaxin Yu
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


    def Squirt_anisotropic(Sdry, Sdry_hp):
        """Predict wet unrelaxed frame compliances at very high frequency from dry frame compliances for transversely isotropic rocks using theoretical formula derived by Mukerji and Mavko, (1994)

        Args:
            Sdry (array):dry rock compliances [S11 S12 S13 S33 S44]
            Sdry_hp (array):dry rock compliances at very high effective stress [S11 S12 S13 S33 S44]

        Returns:
            S_wet (array) : The wet-frame compliances [S11 S12 S13 S33 S44]
        Refs:
        Mukerji and Mavko (1994) Geophysics.
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

    def White_Dutta_Ode(Kdry, Gdry, K0, phi, rho0, rhofl1,rhofl2, Kfl1, Kfl2,eta1,eta2,kapa,a,sg,freq ):
        """Dispersion and Attenuation of partial saturation using White and Dutta–Odé Model. 

        Args:
            Kdry (GPa): bulk modulus of the dry rock 
            Gdry (GPa): shear modulus of the dry rock 
            K0 (GPa): Isotropic mineral bulk modulus
            phi (frac): porosity
            rho0 (g/cm3): mineral density
            rhofl1 (g/cm3): density of the fluid opcupying the central sphere
            rhofl2 (g/cm3): density of the fluid opcupying the outer sphere
            Kfl1 (GPa): bulk modulus of the fluid opcupying the central sphere
            Kfl2 (GPa): bulk modulus of the fluid opcupying the outer sphere
            eta1 (poise): viscousity of the fluid opcupying the central sphere
            eta2 (poise): viscousity of the fluid opcupying the outer sphere
            kapa (cm2): absolute permeability of the rock
            a (cm): radius of central sphere , sg=a3/b3
            sg (frac): saturation of fluid opcupying the central sphere
            freq (scalar or array): frequencies 

            
            
            
        Returns:
            Vp (m/s):
            a_w: 
            K_star: _description_
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
