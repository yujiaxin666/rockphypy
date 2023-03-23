#!/usr/bin/env python
# -*-coding:utf-8 -*-



## Effective contact-based models blended with elasic bounds for calculating the elastic properties of granular medium.

import numpy as np
from rockphypy.utils import utils


class GM: 
    """
    Contact based granular medium models and extensions.
    
    """    
    def ThomasStieber(phi_sand, phi_sh, vsh):
        """Thomas-Stieber porosity model for sand-shale system. 

        Parameters
        ----------
        phi_sand : float
            clean sand porosity
        phi_sh : float
            shale porosity
        vsh : float or array-like
            volume faction of shale in the mixture

        Returns
        -------
        float or array-like
            phi_ABC,phi_AC (frac): porosity line as shown in Fig 5.3.2 in (Mavko,2020)
        """    
    
        phi1=phi_sand-(1-phi_sh)*vsh# dirty sand
        phi2=phi_sh*vsh
        phi_ABC= np.maximum(phi1,phi2)# take values sitting above 
        phi_AD=phi_sand+ phi_sh*vsh # structured shale 
        ##join point A and C
        m= phi_sh-phi_sand
        b= phi_sand
        phi_AC= m*vsh+b
        return phi_ABC,phi_AC  
    def silty_shale(C, Kq,Gq, Ksh, Gsh):
        """Dvorkin–Gutierrez silty shale model: model the elastic moduli of decreasing clay content for shale. 

        Parameters
        ----------
        C : float or array-like
            volume fraction of clay
        Kq : float
            bulk modulus of silt grains
        Gq : float
            shear modulus of silt grains
        Ksh : float
            saturated bulk modulus of pure shale
        Gsh : float
            saturated shear modulus of pure shale, * Ksh and Gsh could be derived from well-log measurements of VP, VS and density in a pure shale zone.

        Returns
        -------
        float or array-like
            K_sat, G_sat: elastic moduli of the saturated silty shale.
        """       

        K_sat = ( C/(Ksh+4*Gsh/3)+ (1-C)/(Kq+4*Gq/3) )**-1 -4*Gsh/3
        Zsh = Gsh/6 *(9*Ksh+8*Gsh)/(Ksh+2*Gsh)
        G_sat = (C/(Gsh+Zsh) + (1-C)/(Gq+Zsh))**-1 -Zsh
        return K_sat,G_sat
    
    def shaly_sand(phis, C, Kss,Gss, Kcc, Gcc):
        """Modeling elastic moduli for sand with increasing clay content using LHS bound rather than using Gassmann relation. 

        Parameters
        ----------
        phis : float
            critical porosity of sand composite
        C : float or array-like
            clay content 
        Kss : float
            saturated bulk moduli for clean sandstone using e.g. HM
        Gss : float
            saturated shear moduli for clean sandstone using e.g. HM
        Kcc : float
            saturated bulk moduli calculated from the sandy shale model at critical clay content using silty shale model
        Gcc : float
            saturated shear moduli calculated from the sandy shale model at critical clay content using silty shale model

        Returns
        -------
        float or array-like
            K_sat,G_sat: saturated rock moduli of the shaly sand 
        """    
    
        K_sat = ( (1-C/phis)/(Kss+4*Gss/3)+ (C/phis)/(Kcc+4*Gss/3) )**-1 -4*Gss/3
        Zss = Gss/6 *(9*Kss+8*Gss)/(Kss+2*Gss)
        G_sat = ((1-C/phis)/(Gss+Zss) + (C/phis)/(Gcc+Zss))**-1 -Zss
        return K_sat,G_sat
    def contactcement(K0, G0, Kc, Gc, phi, phic, Cn,  scheme):
        """Compute dry elastic moduli of cemented sandstone via Contact cement model by Dvorkin &Nur (1996).

        Parameters
        ----------
        K0 : float
            Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
        Kc : float
            Bulk modulus of cement
        Gc : float
            Shear modulus of cement
        phi : float or array-like
            Porosity
        phic : float
            Critical Porosity
        Cn : float
            coordination number
        scheme : int
            Scheme of cement deposition
                    1=cement deposited at grain contacts
                    2=cement deposited at grain surfaces

        Returns
        -------
        _type_
            K_dry, G_dry (GPa): Effective elastic moduli of dry rock

        References
        ----------
        - Dvorkin & Nur, 1996, Geophysics, 61, 1363-1370
        
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

    def hertzmindlin( K0, G0, phic, Cn, sigma, f):
        """Compute effective dry elastic moduli of granular packing under hydrostatic pressure condition via Hertz-Mindlin approach. Reduced shear factor that honours the non-uniform contacts in the granular media is implemented.

        Parameters
        ----------
        K0 : float
            Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
	    phic : float
            Critical Porosity
        Cn : float
            coordination number
        sigma : float or array-like
            effective stress
        f : float
            reduced shear factor between 0 and 1
            0=dry pack with inifinitely rough spheres; 
            1=dry pack with infinitely smooth spheres

        Returns
        -------
        K_dry, G_dry : float or array-like
            effective elastic moduli of dry pack

        References
	    ----------
        - Rock physics handbook section 5.5.
        - Bachrach, R. and Avseth, P. (2008) Geophysics, 73(6), E197–E209.

        """        
    
        sigma =sigma/1000 # converts pressure unit to GPa
        nu=(3*K0-2*G0)/(6*K0+2*G0) # poisson's ratio of mineral mixture
        K_dry = (sigma*(Cn**2*(1-phic)**2*G0**2) / (18*np.pi**2*(1-nu)**2))**(1/3)
        G_dry = ((2+3*f-nu*(1+3*f))/(5*(2-nu))) * ((sigma*(3*Cn**2*(1-phic)**2*G0**2)/(2*np.pi**2*(1-nu)**2)))**(1/3)
        return K_dry, G_dry


    def softsand(K0, G0, phi, phic, Cn, sigma, f):
        """Soft-sand (unconsolidated sand) model: model the porosity-sorting effects using the lower Hashin-Shtrikman-Walpole bound. (Also referred to as the 'friable-sand model' in Avseth et al. (2010). 

        Parameters
        ----------
        K0 : float
            Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
        phi : float or array like
            Porosity
	    phic : float
            Critical Porosity
        Cn : float
            coordination number
        sigma : float or array-like
            effective stress
        f : float
            reduced shear factor between 0 and 1
            0=dry pack with inifinitely rough spheres; 
            1=dry pack with infinitely smooth spheres

        Returns
        -------
        float or array-like
            K_dry, G_dry (GPa): Effective elastic moduli of dry pack
        
        References
	    ----------
        - The Uncemented (Soft) Sand Model in Rock physics handbook section 5.5
        - Avseth, P.; Mukerji, T. & Mavko, G. Cambridge university press, 2010
        """        
           
        K_HM, G_HM = GM.hertzmindlin(K0, G0, phic, Cn, sigma, f)
        K_dry =-4/3*G_HM + (((phi/phic)/(K_HM+4/3*G_HM)) + ((1-phi/phic)/(K0+4/3*G_HM)))**-1
        aux = G_HM/6*((9*K_HM+8*G_HM) / (K_HM+2*G_HM)) # auxiliary variable 
        G_dry = -aux + ((phi/phic)/(G_HM+aux) + ((1-phi/phic)/(G0+aux)))**-1
        return K_dry, G_dry

    def Walton(K0, G0, phic, Cn, sigma, f):
        """Compute dry rock elastic moduli of sphere packs based on the Walton (1987)' thoery. Reduced shear factor that honours the non-uniform contacts in the granular media is implemented. 

        Parameters
        ----------
        K0 : float
            Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
	    phic : float
            Critical Porosity
        Cn : float
            coordination number
        sigma : float or array-like
            effective stress
        f : float
            reduced shear factor between 0 and 1
            0=dry pack with inifinitely rough spheres; 
            1=dry pack with infinitely smooth spheres

        Returns
        -------
        float or array-like
            K_w, G_w: Effective elastic moduli of dry pack
        
        References
	    ----------
        - Walton model in Rock physics handbook section 5.5
        - Walton, K., 1987, J. Mech. Phys. Solids, vol.35, p213-226.
        - Bachrach, R. and Avseth, P. (2008) Geophysics, 73(6), E197–E209
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

    def johnson(K0, G0,n, phi, epsilon, epsilon_axial, path='together'):
        """effective theory for stress-induced anisotropy in sphere packs. The transversely isotropic strain is considered as a combination of hydrostatic strain and uniaxial strain.

        Parameters
        ----------
        K0 : float
            Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
        n : float
            coordination number
        phi : float or array like
            porosity
        epsilon : float or array like
            hydrostatic strain (negative in compression)
        epsilon_axial : float or array like
            uniaxial strain (along 3-axis)
        path : str, optional
            'together': the hydrostatic and uniaxial strains are applied simultaneously
            'uni_iso': the uniaxial strain is applied first followed by a hydrostatic strain
            'iso_uni': the hydrostatic strain is applied first followed by a uniaxial strain by default 'together'

        Returns
        -------
        array and float 
            C: (matrix): VTI stiffness matrix
            sigma33: non zero stress tensor component
            sigma11: non zero stress tensor component, sigma11=sigma22
        
        References
	    ----------
        - Norris, A. N., and Johnson, D. L., 1997, ASME Journal of Applied Mechanics, 64, 39-49.
        - Johnson, D.L., Schwartz, L.M., Elata, D., et al., 1998. Transactions ASME, 65, 380–388.
        
        """        
                
        # factors
        nu0 =  utils.poi(K0, G0)
        Cn = 4*G0/(1-nu0)
        Ct =  (8* G0)/(2-nu0)
        Bw = 2/(np.pi*Cn)
        Cw = 4/np.pi * (1/Ct -1/Cn)
        gamma = 3/32 *Cn*Ct*n*(1-phi)*np.sqrt(-epsilon)
        alpha = np.sqrt(epsilon/epsilon_axial)
        I0 = (1/2)*(np.sqrt(1+alpha**2)+ alpha**2*np.log(1+np.sqrt(1+alpha**2)/alpha))
        I2 = (1/4)*((1+alpha**2)**(3/2)-alpha**2*I0)
        I4 = (1/6)*((1+alpha**2)**(3/2)-3*alpha**2*I2)
        # stiffnesses
        C11 =  gamma/alpha*(4*Bw/3 + 2*Cw/5 + epsilon/epsilon_axial *(2*Bw/15 + Cw/35))
        C33 = (gamma/alpha)* (4*Bw*I2 + 2*Cw*I4)

        C13 = (gamma/alpha)*Cw*(I2-I4)
        C44 = (gamma/alpha)* ((Bw/2)*(I0+I2) + Cw*(I2-I4))
        C66 = (gamma/alpha)* (Bw*(I0-I2) + (Cw/4)*(I0-2*I2+I4))
        # stiffness matrix
        C= utils.write_VTI_matrix(C11,C33,C13,C44,C66)
        
        # the stress component is path dependent
        Temp = -(-epsilon)**(3/2)*(1-phi)*n/ (4*np.pi*alpha**3)
        if path == 'together':# isotropic and uniaxiala strains are applied simultaneously
            sigma33 = 2*Temp*(Ct*(I2-I4)+Cn*(alpha**2*I2+I4))
            sigma11 =Temp*(-Ct*(I2-I4)+Cn*(alpha**2*I0+(1-alpha**2)*I2-I4))  
        elif path == 'uni_iso':
            sigma33 = 2*Temp*(1/12*Ct+Cn*(alpha**2*I2+I4))
            sigma11 =Temp*(-1/12*Ct+Cn*(alpha**2*I0+(1-alpha**2)*I2-I4))
        elif path == 'iso-uni':
            sigma33 = 2*Temp*(Ct*(alpha**2*I0+(1-alpha**2)*I2-I4-2*alpha**3/3)+Cn*(alpha**2*I2+I4))
            sigma11 = Temp*(-Ct*(alpha**2*I0+(1-alpha**2)*I2-I4-2*alpha**3/3+Cn*(alpha**2*I0+(1-alpha**2)*I2-I4)))

        return C, sigma33,sigma11

    # ---------------------more functions added------------------------------#

    def stiffsand(K0, G0, phi, phic, Cn, sigma, f):
        """Stiff-sand model:  Modified Hashin-Shtrikman upper bound with Hertz-Mindlin end point, counterpart to soft sand model. 
        model the porosity-sorting effects using the lower Hashin–Shtrikman–Walpole bound. 

        Parameters
        ----------
        K0 : float
            Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
        phi : float or array like
            Porosity
	    phic : float
            Critical Porosity
        Cn : float
            coordination number
        sigma : float or array-like
            effective stress
        f : float
            reduced shear factor between 0 and 1
            0=dry pack with inifinitely rough spheres; 
            1=dry pack with infinitely smooth spheres

        Returns
        -------
        float or array-like
            K_dry, G_dry (GPa): Effective elastic moduli of dry pack
        """        
          
        K_HM, G_HM = GM.hertzmindlin(K0, G0, phic, Cn, sigma, f)
        K_dry = -4/3*G0 + (((phi/phic)/(K_HM+4/3*G0)) + ((1-phi/phic)/(K0+4/3*G0)))**-1
        aux = G0/6*((9*K0+8*G0) / (K0+2*G0))
        G_dry = -aux + ((phi/phic)/(G_HM+aux) + ((1-phi/phic)/(G0+aux)))**-1
        return K_dry, G_dry

    def constantcement(phi_b, K0, G0, Kc, Gc,phi, phic, Cn,scheme):
        """Constant cement (constant depth) model according to Avseth (2000)

        Parameters
        ----------
        phi_b : _type_
            adjusted high porosity end memeber
        K0 : float
            Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
        Kc : float
            Bulk modulus of cement
        Gc : float
            Shear modulus of cement
        phi : float or array-like
            Porosity
        phic : float
            Critical Porosity
        Cn : float
            coordination number
        scheme : int
            Scheme of cement deposition
                    1=cement deposited at grain contacts
                    2=cement deposited at grain surfaces

        Returns
        -------
        float or array-like
            K_dry, G_dry (GPa): Effective elastic moduli of dry rock

        References
	    ----------
        - Avseth, P.; Dvorkin, J.; Mavko, G. & Rykkje, J. Geophysical Research Letters, Wiley Online Library, 2000, 27, 2761-2764

        """        
            
        K_b,G_b=GM.contactcement(K0, G0, Kc, Gc, phi_b, phic, Cn,  scheme)
        T=phi/phi_b
        Z=(G_b/6)*(9*K_b+8*G_b)/(K_b+2*G_b)
        Kdry=(T/(K_b+4*G_b/3)+(1-T)/(K0+4*G_b/3))**(-1)-4*G_b/3
        Gdry=(T/(G_b+Z)+(1-T)/(G0+Z))**(-1)-Z

        return Kdry, Gdry  
        
    def MUHS(K0, G0, Kc,Gc,phi, phi_b,phic, Cn,scheme):
        """Increasing cement model: Modified Hashin-Strikmann upper bound blend with contact cement model. For elastically stiff sandstone modelling.

        Parameters
        ----------
        K0 : float
            Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
        Kc : float
            Bulk modulus of cement
        Gc : float
            Shear modulus of cement
        phi : float or array-like
            Porosity
        phi_b : _type_
            adjusted high porosity end memeber
        phic : float
            Critical Porosity
        Cn : float
            coordination number
        scheme : int
            Scheme of cement deposition
                    1=cement deposited at grain contacts
                    2=cement deposited at grain surfaces

        Returns
        -------
        float or array-like
            K_dry, G_dry (GPa): Effective elastic moduli of dry rock

        References
	    ----------
        - Avseth, P.; Mukerji, T. & Mavko, G. Cambridge university press, 2010
        
        """        
         
        # adjusted high porosity end memeber 
        K_b, G_b=GM.contactcement(K0, G0, Kc, Gc, phi_b, phic, Cn, scheme)
        # interpolation between the high-porosity end member Kb, Gb and the mineral point. 
        K_dry = -4/3*G0 + (((phi/phi_b)/(K_b+4/3*G0)) + ((1-phi/phi_b)/(K0+4/3*G0)))**-1
        tmp = G0/6*((9*K0+8*G0) / (K0+2*G0))
        G_dry = -tmp + ((phi/phi_b)/(G_b+tmp) + ((1-phi/phi_b)/(G0+tmp)))**-1
        return K_dry, G_dry

    def Digby(K0, G0, phi, Cn, sigma, a_R):
        """Compute Keff and Geff using Digby's model

        Parameters
        ----------
        K0 : float
            Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
        phi : float
            Porosity
        Cn : float
            coordination number
        sigma : float or array-like
            stress
        a_R : float
            a_R (unitless): ratio of the radius of the initially bonded area to the grain radius

        Returns
        -------
        float or array-like
            Keff, Geff (Gpa): effective medium stiffness

        References
	    ----------
        - Digby, P.J., 1981. Journal of Applied Mechanics, 48, 803–808.
        """        
           
        sigma =sigma/1000 # converts pressure unit to GPa
        nu=(3*K0-2*G0)/(6*K0+2*G0) # poisson's ratio of mineral mixture
        
        coeff= np.array([1,0, 3/2* a_R**2, - 3*np.pi*(1-nu)*sigma/(2*Cn*(1-phi)*G0)]) # coeffs of cubic equation
        roots= np.roots(coeff)
        root= roots.real[roots.real>=0]
        b_R= np.sqrt(root**2+a_R**2)
        Sn_R= 4*G0*b_R/(1-nu)
        St_R= 8*G0*a_R/(2-nu)
        Keff= Cn*(1-phi)*Sn_R/12/np.pi
        Geff= Cn*(1-phi)*(Sn_R+1.5*St_R)/20/np.pi
        return Keff, Geff


    def pcm(f,sigma, K0,G0,phi, phic, v_cem,v_ci, Kc,Gc, Cn, mode,scheme,f_):
        """Computes effective elastic moduli of patchy cemented sandstone according to Avseth (2016). 

        Parameters
        ----------
        f : float 
            volume fraction of cemented rock in the binary mixture    
        sigma : float or array-like
            effective stress
        K0 : float
            Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
        phi : float
            Porosity
        phic : float
            Critical Porosity
        v_cem : float
            cement fraction in contact cement model. phi_cem= phic-vcem 
        v_ci : float
            cement threshold above which increasing cement model is applied 
        Kc : float
            bulk modulus of cement
        Gc : float
            shear modulus of cement
        Cn : float
            coordination number
        mode : str
            'stiff' or 'soft'. stiffest mixing or softest mixing. Defaults to 'stiff'.
        scheme : int
            contact cement scheme. 
            1=cement deposited at grain contacts
            2=cement deposited at grain surfaces
        f_ : float
            slip factor in HM modelling. Defaults to 0.5.

        Note
        ----
            (Avseth,2016): If 10% is chosen as the “critical” cement limit, the  increasing cement model can be used in addition to the contact cement model. (Torset, 2020): with the increasing cement model appended at 4% cement"

        Returns
        -------
        float or array-like
            K_DRY, G_DRY (GPa): effective elastic moduli of the dry rock

        References
	    ----------
            - Avseth, P.; Skjei, N. & Mavko, G. The Leading Edge, GeoScienceWorld, 2016, 35, 868-87.

        """        
          
        #------------------------------------------------
        # compute the unconsolidated end member
        #------------------------------------------------
        
        Kunc,Gunc= GM.hertzmindlin(K0, G0, phic=phic, Cn=Cn, sigma=sigma, f=f_)
        #------------------------------------------------
        # compute the cemented end member
        #------------------------------------------------
        if v_cem<=v_ci:
            Kcem,Gcem= GM.contactcement(K0, G0,Kc, Gc, phic-v_cem, phic=phic, Cn=Cn, scheme=scheme)

        #------------------------------------------------
        # increasing cement model 
        #------------------------------------------------
        else: 
            Kcem,Gcem= GM.MUHS(K0, G0, Kc,Gc,phic-v_cem, phic-v_ci,phic=phic, Cn=Cn, scheme=scheme)

        #------------------------------------------------
        # first HS 
        #------------------------------------------------   

        if mode == 'stiff': # stiffest mixing,shelf is cemented sand, inner core is unconsolidated sand  
            Kp=Kcem+ (1-f)/( (Kunc-Kcem)**-1 + f*(Kcem+4*Gcem/3)**-1 )

            Temp= (Kcem+2*Gcem)/(5*Gcem *(Kcem+4*Gcem/3))
            Gp=Gcem+(1-f)/( (Gunc-Gcem)**-1 + 2*f*Temp)
        else:  # soft mixing inner core is cemented sand, outer shelf is unconsolidated sand  
            Kp=Kunc+ f/( (Kcem-Kunc)**-1 + (1-f)*(Kunc+4*Gunc/3)**-1 )

            Temp= (Kunc+2*Gunc)/(5*Gunc *(Kunc+4*Gunc/3))
            Gp=Gunc+f/( (Gcem-Gunc)**-1 + 2*(1-f)*Temp)
        
        #------------------------------------------------Second HS: interpolate the effective high-porosity end member and the mineral point using MLHS bound 
        #------------------------------------------------
        K_DRY =-4/3*Gp + (((phi /phic)/(Kp+4/3*Gp)) + ((1-phi/phic)/(K0+4/3*Gp)))**-1
        tmp = Gp/6*((9*Kp+8*Gp) / (Kp+2*Gp))
        G_DRY = -tmp + ((phi/phic)/(Gp+tmp) + ((1-phi/phic)/(G0+tmp)))**-1

        return K_DRY, G_DRY

