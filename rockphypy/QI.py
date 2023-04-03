import numpy as np
from rockphypy.utils import utils
from rockphypy.EM import EM
from rockphypy.GM import GM
from rockphypy.Fluid import Fluid

class QI:
    """
    Useful functionalities for quantitative intepretation and relevant tasks.
    """ 

    @staticmethod 
    def matrix_modulus(vsh,phi_c,phi_0,M_sh,M_g, M_c):
        """Calculate the modulus of rock matrix as a function of porosity variation caused by cementation, Note that the rock matrix contains everything excluding pore space.

        Parameters
        ----------
        vsh : float 
            bulk volume fraction of shale, can be derived from well log.
        phi_c : float
            critical porosity
        phi_0 : float or array
            static porosity during cementation ranging between 0.4 to 0 should be phi when phi is an array of porosity
        M_sh : float
            modulus of shale
        M_g : float
            modulus of grain material
        M_c : float
            modulus of cement

        Returns
        -------
        float or array
            M_mat: updated matrix modulus 
        """    

        #phi_0=1-np.linspace(0.0,phi_c,100) # porosity state from the begining of cementation to the end of cementation 
        fc = phi_c - phi_0 # Fraction of cement in the rock
        
        vsh_s=vsh/(1-phi_0)# shale volume fraction in matrix
        fgs = (1-phi_c-vsh)/(1-phi_0) #faction of grain in the matrix
        fcs = fc/(1-phi_0) # fraction of cement  in the matrix
        if vsh==0:
            volume=np.vstack((fgs,fcs)).T
            M=np.array([M_g,M_c])
            _,_,M_mat=EM.VRH(volume, M)
            return M_mat
        else:
            volume=np.vstack((vsh_s,fgs,fcs)).T
            M=np.array([M_sh, M_g,M_c])
            _,_,M_mat=EM.VRH(volume, M)  

            return M_mat ## it is a numpy ndarray

    @staticmethod 
    def den_matrix(vsh,phi_c,phi_0,den_sh,den_g, den_c):
        """Calculate the matrix density as a function of porosity variation caused by cementation.

        Parameters
        ----------
        vsh : float 
            bulk volume fraction of shale, can be derived from well log.
        phi_c : float
            critical porosity
        phi_0 : float or array
            static porosity during cementation ranging between 0.4 to 0 should be phi when phi is an array of porosity
        den_sh : float
            density of the clay 
        den_g : float
            density of the grain
        den_c : float
            density of the cement

        Returns
        -------
        float or array
            den_mat: updated matrix density
        """    

        fc = phi_c - phi_0 # Fraction of cement in the rock
        
        vsh_s=vsh/(1-phi_0)# shale volume fraction in matrix
        fgs = (1-phi_c-vsh)/(1-phi_0) #faction of grain in the matrix
        fcs = fc/(1-phi_0) # fraction of cement  in the matrix
        volume=np.vstack((vsh_s,fgs,fcs)).T
        D=np.array([den_sh,den_g, den_c])
        den_mat,_,_=EM.VRH(volume, D)
        
        return den_mat

    @staticmethod 
    def screening(Dqz,Kqz,Gqz,Dsh,Ksh,Gsh,Dc,Kc,Gc,Db,Kb,phib_p,phi_c,sigma,vsh,scheme,f, Cn):
        """compute elastic bounds used for rock physics screening, the lower bound is computed using friable sand model, and the upper bound is contact cement model blend with increasing cement model.  

        Parameters
        ----------
        K0 : float
                Bulk modulus of grain material in GPa
        G0 : float
            Shear modulus of grain material in GPa
        Dqz : float
            Density of the grain. not limited to quartz grain
        Kqz : float
            Bulk modulus of grain material in GPa
        Gqz : float
            Shear modulus of grain material in GPa
        Dsh : float
            density the clay
        Ksh : float
            bulk modulus of the clay
        Gsh : float
            shear modulus of the clay
        Dc : float
            density of the cement 
        Kc : float
            Bulk modulus of cement
        Gc : float
            Shear modulus of cement
        Db : float
            density of the pore fluid
        Kb : float
            bulk modulus of the pore fluid         
        phib_p : float
            adjusted high porosity end memeber
        phic : float
                Critical Porosity
        sigma : float or array-like
            effective stress
        vsh : _type_
            _description_
        scheme : int
            Scheme of cement deposition
                    1=cement deposited at grain contacts
                    2=cement deposited at grain surfaces
        f : float
            reduced shear factor between 0 and 1
                0=dry pack with inifinitely rough spheres; 
                1=dry pack with infinitely smooth spheres
        Cn : float
            coordination number

        Returns
        -------
        array
            phi,vp1,vp2,vp3,vs1,vs2,vs3: porosity and velocities required for elastic diagnostics bounds
        """    
        
        #------------calculate the effective bulk and shear modulus------------#
        _,_,K0=EM.VRH(np.array([vsh,1-vsh]),np.array([Ksh,Kqz]))
        _,_,G0=EM.VRH(np.array([vsh,1-vsh]),np.array([Gsh,Gqz])) # this is only two component sandstone, when sandstone is multimineralogy, vectorization should apply here
        D0=Dsh*vsh+Dqz*(1-vsh)
        
        phi = np.linspace(1e-7,phi_c,100) #define porosity range according to critical porosity
    
        Kdry1, Gdry1 = GM.softsand(K0, G0, phi, phi_c, Cn, sigma,f) # frible sand 
    
        Kdry3, Gdry3= GM.contactcement(K0, G0, Kc, Gc, phi, phi_c, Cn, scheme) # contact cement

        K_mat=QI.matrix_modulus(vsh,phi_c,phi,Ksh,Kqz,Kc) # update the mineral matrix modulus for contact cement model
        Den_mat=QI.den_matrix(vsh,phi_c,phi,Dsh,Dqz,Dc) # update matrix density for contact cement model
        
        Ksat1, Gsat1 = Fluid.Gassmann(Kdry1,Gdry1,K0,Kb,phi) 
    
        Ksat3, Gsat3 = Fluid.Gassmann(Kdry3,Gdry3,K_mat,Kb,phi) 

        rho  = D0*(1-phi)+Db*phi # update bulk density for saturated soft sand model and stiff sand
        rho_cst=Den_mat*(1-phi)+Db*phi # update bulk density for saturated contact cement model
        vp1,vs1= utils.V(Ksat1,Gsat1,rho)

        vp3,vs3= utils.V(Ksat3,Gsat3,rho_cst)
        
        vp3[phi<phib_p]=np.nan # contact cement model is invalid for small porosity, here we set the critical cement limit is 10$
        vs3[phi<phib_p]=np.nan

        Kdry2, Gdry2 = GM.MUHS(K0, G0, Kc,Gc,phi, phib_p,phi_c, Cn,scheme) #increasing cement model 
        Ksat2, Gsat2 = Fluid.Gassmann(Kdry2,Gdry2,K0,Kb,phi) 
        vp2,vs2= utils.V(Ksat2,Gsat2,rho)
        vp2[phi>phib_p]=np.nan # increasing cement model is invalid for small porosity
        vs2[phi>phib_p]=np.nan

        return phi,vp1,vp2,vp3,vs1,vs2,vs3
