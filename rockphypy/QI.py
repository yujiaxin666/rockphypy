
import matplotlib.colors
cmap_sand = matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","orange","yellow"])
cmap_kde=matplotlib.colors.LinearSegmentedColormap.from_list("", ['white','grey','#4895d9','#65d698',"yellow",'yellow'])
# colormap  
cmap_cem=matplotlib.colors.LinearSegmentedColormap.from_list("", [(0,'blue'),(0.05,'#0a75ad'),(0.05,'#4ca3dd'),(0.1,'#20b2aa'),(0.4,"#008000"),(0.6,"yellow"),(1,"red")])

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from kde_diffusion import kde2d

# import matplotlib.pyplot as plt


from rockphypy.utils import utils
from rockphypy.EM import EM
from rockphypy.GM import GM
from rockphypy.Emp import Empirical
from rockphypy.Fluid import Fluid


class QI:
    
    """
    Useful functionalities for quantitative intepretation and relevant tasks.
    """ 
    def __init__(self,Vp,phi,**kwargs): #phi,den,Vsh,eff_stress, TVD,
        """initialize the parameters for various QI plots. e.g., phi,den,Vsh,eff_stress, TVD can be given depending on required the input parameters to the plot funtion.

        Parameters
        ----------
        Vp : array
            p wave velocity 
        phi : array
            porosity 
        """        
        self.Vp=Vp
        self.phi=phi
        self.Vsh=kwargs.get('Vsh')
        self.Vs=kwargs.get('Vs')
        self.den=kwargs.get('den')
        self.eff_stress= kwargs.get('eff_stress')
        self.TVD=kwargs.get('TVD')
        pass
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
    def screening_plot(self,phi,vp1,vp2,vp3,cmap=cmap_sand):
        """plot the rock physics screening crossplot 

        Parameters
        ----------
        phi : array
            porosity 
        vp1 : array
            lower bound velocities modeled using friable sandstone model by default . 
        vp2 : array
            upper bound velocities modeled using MUSH model by default. 
        vp3 : array
            upper bound velocities modeled using contact cement model by default. 
        cmap : string, optional
            colormap, can be default colormaps in matplotlib, by default using the customized colormap: cmap_sand
        Returns
        -------
        object
            elastic bounds plot used for rock physics screening 
        """      

        fig,ax=plt.subplots()
        fig.set_size_inches(7, 6)
        ax.plot(phi,vp3,'-k', lw=4, alpha=0.7)
        ax.plot(phi,vp1,'--k', lw=2, alpha=0.7)
        ax.plot(phi,vp2,'-k',lw=4, alpha=0.7)
        ax.set_ylabel('Vp')
        ax.set_xlabel('Porosity')
        ax.grid(ls='--',alpha=0.7)            
        plt.scatter(self.phi,self.Vp,c=1-self.Vsh,vmin=0, vmax=1,edgecolors='grey',s=100,alpha=1,cmap=cmap)
        cbar=plt.colorbar()
        cbar.set_label(r'$\rm V_{Sand}$')
        return fig
    
    @staticmethod
    def normalize(density):

        '''
        normalize the kde with respect to the maximum value and map value to 0-1
        '''

        den_max=density.max()
        
        # map value
        new=1-(den_max-density)/den_max
        return new

    def kde_plot(self,phi,v1,v2,v3,cmap=cmap_kde,vels='Vp',n=64):
        """plot field data or measurements as 2D probability density functions in the elastic bounds cross plot

        Parameters
        ----------
        phi : array
            porosity 
        vp1 : array
            lower bound velocities modeled using friable sandstone model by default . 
        vp2 : array
            upper bound velocities modeled using MUSH model by default. 
        vp3 : array
            upper bound velocities modeled using contact cement model by default. 
        cmap : string, optional
            colormap, can be default colormaps in matplotlib, by default using the customized colormap: cmap_kde
        vels : str, optional
            choose either P wave or S wave velocity for plotting, by default 'Vp'
        n : int, optional
            grid parameter used in KDE-diffusion, by default 64

        Returns
        -------
        object
            KDE plot with elastic bounds
        """        
        #KDE via diffusion on data
        fig,ax=plt.subplots()
        fig.set_size_inches(7, 6)
        ax.plot(phi,v3,'-k', lw=4, alpha=0.7)
        ax.plot(phi,v1,'--k', lw=2, alpha=0.7)
        ax.plot(phi,v2,'-k',lw=4, alpha=0.7)
        #ax.set_ylabel('Vp')
        ax.set_xlabel('Porosity')
        
        ax.set_xlim([-0.007,0.407])
        if vels=='Vp':
            (density, grid, bandwidth) = kde2d(self.phi, self.Vp, n=n, limits=None)
            # map density to 0 and 1
            density=QI.normalize(density)
            # extremely low density is mapped to 0
            density[density<0.0001]=0
            plt.pcolormesh(grid[0],grid[1], density.T, cmap=cmap_kde,shading='auto')
            cbar=plt.colorbar()
            cbar.set_label('KDE Normalized')
            plt.contour(grid[0],grid[1], density.T,colors='white',linewidths=0.5)
            ax.set_ylabel('Vp')
            
        elif vels=='Vs':
            (density, grid, bandwidth) = kde2d(self.phi, self.Vs, n=n, limits=None)
            # map density to 0 and 1
            density=QI.normalize(density)
            # extremely low density is mapped to 0
            density[density<0.0001]=0
            plt.pcolormesh(grid[0],grid[1], density.T, cmap=cmap_kde,shading='auto')
            cbar=plt.colorbar()
            cbar.set_label('KDE Normalized')
            plt.contour(grid[0],grid[1], density.T,colors='white',linewidths=0.5)
            ax.set_ylabel('Vs')
        ax.grid(alpha=1)
        return fig
    @staticmethod
    def cst_vels(phi_b,K0,D0,G0,phi,phi_c,Cn,Kc,Gc,Db,Kb,scheme,vsh,Dsh,Dqz,Dc):
        """compute velocities using constant cement model at different cement amounts 

        Parameters
        ----------
        phi_b : float
            adjusted porosity for constant cemnet model 
        K0 : float
            Bulk modulus of grain material in GPa
        D0 : float
            Density of grain material
        G0 : float
            Shear modulus of grain material in GPa
        phi : float
            porosity
        phi_c : float
            Critical Porosity
        Cn : float
            critical porosity
        Kc : float
            Bulk modulus of cement
        Gc : float
            Shear modulus of cement
        Db : float
            density of the pore fluid
        Kb : float
            bulk modulus of the pore fluid 
        scheme : int
            Scheme of cement deposition
                    1=cement deposited at grain contacts
                    2=cement deposited at grain surfaces
        vsh : float
            shale content
        Dsh : float
            density the clay
        Dqz : float
            Density of the grain. not limited to quartz grain
        Dc : float
            density of the cement

        Returns
        -------
        array
            vp,vs: velocities given by constant cement model 
        """        
        Kdry, Gdry=GM.constantcement(phi_b, K0, G0,Kc, Gc, phi, phi_c, Cn,scheme )
        D=QI.den_matrix(vsh,phi_c,phi_b,Dsh,Dqz,Dc) 
        vp,vs,_= Fluid.vels(Kdry,Gdry,K0,D,Kb,Db,phi)
        vp[phi>phi_b]=np.nan
        vs[phi>phi_b]=np.nan
        return vp,vs

    @staticmethod
    def cst_plot(Dqz,Kqz,Gqz,Dsh,Ksh,Gsh,Dc,Kc,Gc,Db,Kb,phib,phib_p,phi_c,sigma,vsh,Cn, scheme,f):
        """Diagnostic plot with constant cement model lines 

        Parameters
        ----------
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
        phib : float
            adjusted high porosity end memeber for constant cement model 
        phib_p : array or list 
            posoities used to drawing the constant cement lines
        phi_c : float
                Critical Porosity
        sigma : float
            effective stress
        vsh : float
            shale content
        Cn : float
            coordination number
        scheme : int
            Scheme of cement deposition
                    1=cement deposited at grain contacts
                    2=cement deposited at grain surfaces
        f : float
            reduced shear factor between 0 and 1
                0=dry pack with inifinitely rough spheres; 
                1=dry pack with infinitely smooth spheres

        Returns
        -------
        object
            fig,ax: constant cement diagnostic plot
        """        
        _,_,K0=EM.VRH(np.array([vsh,1-vsh]),np.array([Ksh,Kqz]))
        _,_,G0=EM.VRH(np.array([vsh,1-vsh]),np.array([Gsh,Gqz])) # this is only two component sandstone, when sandstone is multimineralogy, vectorization should apply here
        D0=Dsh*vsh+Dqz*(1-vsh)
        
        phi,vp1,vp2,vp3,vs1,vs2,vs3 = QI.screening(Dqz,Kqz,Gqz,Dsh,Ksh,Gsh,Dc,Kc,Gc,Db,Kb,phib,phi_c,sigma,vsh,scheme,f, Cn)

        fig,ax=plt.subplots()
        ax.plot(phi,vp3,'-k', lw=4, alpha=0.7)
        ax.plot(phi,vp1,'-k', lw=3, alpha=0.7)

        for i in phib_p:
            
            ax.plot(phi,QI.cst_vels(i,K0,D0,G0,phi,phi_c,Cn,Kc,Gc,Db,Kb,scheme,vsh,Dsh,Dqz,Dc)[0],'--k',lw=2, alpha=0.7)
        ax.set_ylabel('Vp (m/s)')
        ax.set_xlabel('porosity')
        ax.grid()
        #plt.legend()
        return fig,ax
    @staticmethod
    def cal_v_const(Dqz,Kqz,Gqz,Dsh,Ksh,Gsh,Dc,Kc,Gc,Db,Kb,phi_b,phi_c,vsh,phi,scheme):
        """input real data porosity and caculate the theoretical constant cement model velocity value.  Note: input porosity cannot be zero, otherwise the returned velocities are Nan.

        Parameters
        ----------
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
        phi_b : _type_
            _description_
        phi_c : float
                Critical Porosity
        vsh : float
            shale content
        phi : array
            porosity
        scheme : int
            Scheme of cement deposition
                    1=cement deposited at grain contacts
                    2=cement deposited at grain surfaces

        Returns
        -------
        array
            vp,vs: constant cement velocities 
        """        
       
        _,_,K0=EM.VRH(np.array([vsh,1-vsh]),np.array([Ksh,Kqz]))
        _,_,G0=EM.VRH(np.array([vsh,1-vsh]),np.array([Gsh,Gqz])) # this is only two component sandstone, when sandstone is multimineralogy, vectorization should apply here
        D0=Dsh*vsh+Dqz*(1-vsh)
        Cn=Empirical.Cp(phi_c) ## calculate coordination number 
        Kdry, Gdry=GM.constantcement(phi_b, K0, G0, Kc, Gc, phi, phi_c, Cn, scheme )
        vp,vs,_= Fluid.vels(Kdry,Gdry,K0,D0,Kb,Db,phi)
        return vp,vs 


    def estimate_cem(self,vcem_seeds,Kqz,Gqz,Ksh,Gsh,phi_c,Cn,Kc,Gc,Db,Kb,scheme,vsh,Dsh,Dqz,Dc):
        """Estimate cement amount for well log data using constant cement model crossplot. 

        Parameters
        ----------
        vcem_seeds : array or list
            some predefined values at which constant cement lines are calculated 
        Kqz : float
            Bulk modulus of grain material in GPa
        Gqz : float
            Shear modulus of grain material in GPa
        Ksh : float
            bulk modulus of the clay
        Gsh : float
            shear modulus of the clay
        phi_c : float
                Critical Porosity
        Cn : float
            coordination number
        Kc : float
            Bulk modulus of cement
        Gc : float
            Shear modulus of cement
        Db : float
            density of the pore fluid
        Kb : float
            bulk modulus of the pore fluid        
        scheme : int
            Scheme of cement deposition
                    1=cement deposited at grain contacts
                    2=cement deposited at grain surfaces
        vsh : float
            shale content
        Dsh : float
            density the clay
        Dqz : float
            Density of the grain. not limited to quartz grain
        Dc : float
            density of the cement 

        Returns
        -------
        array
            cement amount estimation for each well log data points
        """        
        _,_,K0=EM.VRH(np.array([vsh,1-vsh]),np.array([Ksh,Kqz]))
        _,_,G0=EM.VRH(np.array([vsh,1-vsh]),np.array([Gsh,Gqz])) # this is only two component sandstone, when sandstone is multimineralogy, vectorization should apply here
        D0=Dsh*vsh+Dqz*(1-vsh)
        # create an empty dataframe 
        df = pd.DataFrame()
        # corresponding porosity is phi_c-vcem_seeds[i]
        for i in vcem_seeds:
            
            VP,VS=QI.cal_v_const(Dqz,Kqz,Gqz,Dsh,Ksh,Gsh,Dc,Kc,Gc,Db,Kb,phi_c-i,phi_c,vsh,self.phi,scheme)
            df['vcem_vp'+str(i)]=VP
        # decide if Vcem is 0 or not, otherwise, interpolate between two closet velocities
        df['vcem']=np.nan
        for i, val in enumerate(self.Vp):
            array= np.array(df.iloc[i, :-1])
            if val <= array[0]:
                df.vcem.values[i]=0
                #df.iloc[i]['vcem']=0
            elif val >=array[-1]:
                df.vcem.values[i]=0.1
                #df.iloc[i]['vcem']=0.1
            else:
                big =  array[array<val].max()# the largest element of myArr less than myNumber
                small = array[array>val].min()# the smallest element of myArr greater than myNumber
                #df.vcem.values[i]
                df.vcem.values[i]=(val-small)/(big-small)*(vcem_seeds[array==big]-vcem_seeds[array==small])+vcem_seeds[array==small]
        return df.vcem.values
    def cement_diag_plot(self,vcem,Dqz,Kqz,Gqz,Dsh,Ksh,Gsh,Dc,Kc,Gc,Db,Kb,phib,phib_p,phi_c,sigma,vsh,Cn, scheme,f):
        """_summary_

        Parameters
        ----------
        vcem : _type_
            _description_
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
        phib : float
            adjusted high porosity end memeber for constant cement model 
        phib_p : array or list 
            posoities used to drawing the constant cement lines
        phi_c : float
                Critical Porosity
        sigma : float or array-like
            effective stress
        vsh : float
            shale content
        Cn : float
            coordination number
        scheme : int
            Scheme of cement deposition
                    1=cement deposited at grain contacts
                    2=cement deposited at grain surfaces
        f : float
            reduced shear factor between 0 and 1
                0=dry pack with inifinitely rough spheres; 
                1=dry pack with infinitely smooth spheres

        Returns
        -------
        object
            cross plot for cement estimation
        """        

        fig,ax=QI.cst_plot(Dqz,Kqz,Gqz,Dsh,Ksh,Gsh,Dc,Kc,Gc,Db,Kb,phib,phib_p,phi_c,sigma,vsh,Cn, scheme,f)
        fig.set_size_inches(7, 6)
        plt.scatter(self.phi,self.Vp,c=vcem,edgecolors='grey',s=100,vmin=0, vmax=0.1,alpha=0.7,cmap=cmap_cem)
        cbar=plt.colorbar()
        cbar.set_label(r'$\rm V_{CEM}$')
        plt.xlabel('Porosity')
        plt.ylim([1500,6500])
        plt.xlim([-0.007,0.507])
        return fig

    # RPT plot 
    @staticmethod
    def plot_rpt(Kdry,Gdry,K0,D0,Kb,Db,Khc,Dhc,phi,sw):
        """Create RPT plot given computed Impedance and Vp/Vs ratio. 

        Parameters
        ----------
        Kdry : float or array
            effective bulk modulus given by rock physics model
        Gdry : float or array
            effective shear modulus given by rock physics model
        K0 : float
            bulk modulus of grain
        D0 : float
            density of grain
        Kb : float
            bulk moduluf of brine 
        Db : float
            density of brine
        Khc : float
            bulk modulus of HC
        Dhc : float
            density of HC
        phi : float or array 
            porosity
        sw : float or array
            water saturation

        Returns
        -------
        python onject: fig
            rpt plot
        """    
        
        # setup empty arrays to store Ip and Vp/Vs values
        IP=np.empty((phi.size,sw.size))
        PS=np.empty((phi.size,sw.size))

        ## loop over Sw, computes elastic moduli of fluid mixture and saturated rock properties with Gassmann's equation
        
        for i,val in enumerate(sw):
            Kf=(val/Kb+(1-val)/Khc)**-1
            #Kf_mix(val,Kb,Khc)
            Df = val*Db+(1-val)*Dhc
            vp,vs,rho= Fluid.vels(Kdry,Gdry,K0,D0,Kf,Df,phi)
            IP[:,i]=vp*rho
            PS[:,i]=vp/vs
        # plot
        fig=plt.figure(figsize=(10,8))
        plt.plot(IP.T, PS.T, '-ok', mec='k', ms=10, mfc='yellow')
        plt.plot(IP[:,-1], PS[:,-1], '-ok', mec='k', ms=10, mfc='blue')# Brine line 

        plt.xlabel('Acoustic Impedance'), plt.ylabel('Vp/Vs')
        # for i,val in enumerate(phi):
        #     plt.text(IP[i,-1],PS[i,-1]+.03,'{:.02f}'.format(val), alpha=1,backgroundcolor='0.9')
        # for i,val in enumerate(sw):
        #     plt.text(IP[-1,i]-100,PS[-1,i],'Gas={:.02f}'.format(1-sw[i]),ha='right',alpha=1)
        return fig