from rockphypy import EM
from rockphypy import Fluid

# specify model parameters
phi=np.linspace(0,1,100,endpoint=True) # solid volume fraction = 1-phi
K0, G0= 37,44 # moduli of grain material
Kw, Gw= 2.2,0 # moduli of water
# VRH bounds
volumes= np.vstack((1-phi,phi)).T
M= np.array([K0,Kw])
K_v,K_r,K_h=EM.VRH(volumes,M)
# Hashin-Strikmann bound
K_UHS,G_UHS= EM.HS(1-phi, K0, Kw,G0,Gw, bound='upper')
# Critical porosity model
phic=0.4 # Critical porosity
phi_=np.linspace(0.001,phic,100,endpoint=True) # solid volume fraction = 1-phi
K_dry, G_dry= EM.cripor(K0, G0, phi_, phic)# Compute dry-rock moduli
Ksat, Gsat = Fluid.Gassmann(K_dry,G_dry,K0,Kw,phi_)# saturate rock with water

# plot
plt.figure(figsize=(6,6))
plt.xlabel('Porosity')
plt.ylabel('Bulk modulus [GPa]')
plt.title('V, R, VRH, HS bounds')
plt.plot(phi, K_v,label='K Voigt')
plt.plot(phi, K_r,label='K Reuss = K HS-')
plt.plot(phi, K_h,label='K VRH')
plt.plot(phi, K_UHS,label='K HS+')
plt.plot(phi_, Ksat,label='K CriPor')
plt.legend(loc='best')
plt.grid(ls='--')