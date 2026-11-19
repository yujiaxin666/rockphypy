#from import_path import rockphypy
from rockphypy import Fluid
import pytest
import numpy as np

class TestFluid:


    @pytest.fixture
    def fluid_input(self):
        vpdry=2.500
        vsdry=1.300
        K0=37
        mu0=46
        rho0=2.75
        rhofl=1.2
        Kfl=2.3
        nu=0.1
        phi=0.2
        perm=0.001
        a=0.02
        alpha=2
        d1=-1
        d2=3
        freq= np.logspace(d1,d2,4,endpoint=True)
        rhodry=(1-phi)*rho0
        Gdry=rhodry*vsdry**2
        Kdry=rhodry*vpdry**2-(4/3)*Gdry
        return Kdry,Gdry,K0,Kfl,rho0,rhofl,nu,phi,perm,a,alpha,freq
    def test_Brie(self):
        Kw, Kgas,e= 2.1,1.2,4
        Sw= np.linspace(0,1,10)
        result = Fluid.Brie(Kw, Kgas, Sw, e)
        expected_result = np.array([1.2       , 1.20013717, 1.20219479, 1.21111111, 1.2351166 , 1.28573388, 1.37777778, 1.52935528, 1.76186557, 2.1       ])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Biot(self,fluid_input):
        Kdry,Gdry,K0,Kfl,rho0,rhofl,nu,phi,perm,a,alpha,freq= fluid_input
        
        result = Fluid.Biot(Kdry,Gdry,K0,Kfl,rho0,rhofl,nu,phi,perm,a,alpha,freq)
        expected_result = (np.array([2.82259127, 2.82512442, 2.82636373, 2.82636977]), np.array([0.28457234, 0.74286946, 0.77719293, 0.77747519]), np.array([1.23458051, 1.25693158, 1.26585154, 1.26589411]), np.array([1.77718961e-04, 1.26089932e-03, 8.71068103e-05, 5.67820847e-06]), np.array([1.39142536e+01, 6.44535023e-01, 2.99632580e-02, 1.95001446e-03]), np.array([3.68817365e-03, 2.27117947e-02, 1.47416907e-03, 9.60675136e-05]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Biot_HF(self,fluid_input):
        Kdry,Gdry,K0,Kfl,rho0,rhofl,_,phi,_,_,alpha,_= fluid_input
        result = Fluid.Biot_HF(Kdry,Gdry,K0,Kfl,rho0,rhofl,phi,alpha)
        expected_result = (2.8263748996477673, 0.7779580777841519, 1.2659329393362635)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Geertsma_Smit_HF(self,fluid_input):
        Kdry,Gdry,K0,Kfl,rho0,rhofl,_,phi,_,_,alpha,_= fluid_input
        result = Fluid.Geertsma_Smit_HF(Kdry,Gdry,K0,Kfl,rho0,rhofl,phi,alpha)
        expected_result =(2.9314866269776054, 1.2659329393362635)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Geertsma_Smit_LF(self,fluid_input):
        Vp0,Vpinf= 2500, 3000
        _,_,_,_,_,rhofl,nu,phi,perm,_,_,freq= fluid_input
        result= Fluid.Geertsma_Smit_LF(Vp0,Vpinf, freq,phi, rhofl, perm, nu)
        expected_result =np.array([2501.12305747, 2754.93967911, 2998.96267556, 2999.99776048])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Gassmann(self,fluid_input):
        Kdry,Gdry,K0,Kfl,_,_,_,phi,_,_,_,_= fluid_input
        result = Fluid.Gassmann(Kdry,Gdry,K0,Kfl,phi)
        expected_result =(14.481969751202783, 3.718000000000001)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Gassmann_sub(self):
        Ksat1= np.arange(15,35,5)
        Kfl1= 2.3
        Kfl2= 1.4
        K0=37
        phi=0.2
        result= Fluid.Gassmann_sub(phi, K0, Ksat1,Kfl1,Kfl2)
        expected_result =np.array([13.0832864 , 18.8777377 , 24.45146108, 29.81683308])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_vels(self,fluid_input):
        Kdry,Gdry,K0,Kfl,_,_,_,phi,_,_,_,_= fluid_input
        D0=2.65
        Dfl= 1.2
        result = Fluid.vels(Kdry,Gdry,K0,D0,Kfl,Dfl,phi)
        expected_result =(2870.0161728517414, 1255.1588460484033, 2.3600000000000003)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Gassmann_vels(self):
        Kfl1= 2.3
        Kfl2= 1.4
        D0=2.65
        Dfl1= 1.2
        Dfl2= 0.8
        Vp1= np.linspace(3500,4000,5)
        Vs1= np.linspace(1500,3000,5)
        rho1= 2.3
        K0=37
        phi=0.2
        result = Fluid.Gassmann_vels(Vp1,Vs1,rho1,Dfl1,Kfl1,Dfl2,Kfl2,K0,phi)
        expected_result =(np.array([1763.02252647, 2203.76025   , 2644.5006265 , 3085.24251897,
        3525.98535888]), np.array([1526.78783106, 1908.48478883, 2290.18174659, 2671.87870436,
        3053.57566213]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Gassmann_approx(self):
        phi=0.2
        Msat1, M0, Mfl1, Mfl2= 35,65,4,3
        result = Fluid.Gassmann_approx(Msat1,M0,Mfl1,phi,Mfl2)
        expected_result =33.76101321585903
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    @pytest.fixture
    def bk_input(self):
        phi=0.2
        Kfl=2.3
        K0, G0= 37,44
        C = np.array([[100, 30, 15,  0,  0,  0],
                    [30, 100, 15,  0,  0,  0],
                    [15, 15, 30,  0,  0,  0],
                    [0,  0,  0, 20,  0,  0],
                    [0,  0,  0,  0, 20,  0],
                    [0,  0,  0,  0,  0, 10]])
        S = np.linalg.inv(C)
        return C,S,K0,G0,Kfl,phi
    def test_Brown_Korringa_dry2sat(self,bk_input):
        _,S,K0,G0,Kfl,phi = bk_input
        result = Fluid.Brown_Korringa_dry2sat(S,K0,G0,Kfl,phi)
        expected_result =np.array([[ 0.01125505, -0.00303066, -0.00333797,  0.        ,  0.        ,
         0.        ],
       [-0.00303066,  0.01125505, -0.00333797,  0.        ,  0.        ,
         0.        ],
       [-0.00333797, -0.00333797,  0.0333532 ,  0.        ,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.05      ,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.05      ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         0.1       ]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Brown_Korringa_sat2dry(self,bk_input):
        _, S,K0,G0,Kfl,phi = bk_input
        result = Fluid.Brown_Korringa_sat2dry(S,K0,G0,Kfl,phi)
        expected_result = np.array([[ 0.01179714, -0.00248858, -0.0056612 ,  0.        ,  0.        ,
         0.        ],
       [-0.00248858,  0.01179714, -0.0056612 ,  0.        ,  0.        ,
         0.        ],
       [-0.0056612 , -0.0056612 ,  0.0433099 ,  0.        ,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.05      ,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.05      ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         0.1       ]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Brown_Korringa_sub(self,bk_input):
        C,_,K0,G0,_,phi = bk_input
        Kfl1,Kfl2= 2.3,1.4
        result = Fluid.Brown_Korringa_sub(C,K0,G0,Kfl1,Kfl2,phi)
        expected_result =(np.array([[99.53530254, 29.53530254, 15.6970462 ,  0.        ,  0.        ,
          0.        ],
        [29.53530254, 99.53530254, 15.6970462 ,  0.        ,  0.        ,
          0.        ],
        [15.6970462 , 15.6970462 , 28.95443071,  0.        ,  0.        ,
          0.        ],
        [ 0.        ,  0.        ,  0.        , 20.        ,  0.        ,
          0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        , 20.        ,
          0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         10.        ]]), np.array([[ 0.01160511, -0.0026806 , -0.00483824,  0.        ,  0.        ,
          0.        ],
        [-0.0026806 ,  0.01160511, -0.00483824,  0.        ,  0.        ,
          0.        ],
        [-0.00483824, -0.00483824,  0.03978293,  0.        ,  0.        ,
          0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.05      ,  0.        ,
          0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.05      ,
          0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
          0.1       ]]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Mavko_Jizba(self):
        Vp_hs, Vs_hs = 4500,3500
        Vpdry, Vsdry = np.arange(2500,4500,500), np.arange(1500,3500,500), 
        K0, rhodry, rhofl,Kfl, phi = 37,2.4,1.2,2.1,0.2
        result = Fluid.Mavko_Jizba(Vp_hs, Vs_hs,Vpdry, Vsdry, K0, rhodry, rhofl,Kfl, phi)
        expected_result = (14.458807159572448, np.array([ 5.57519432,  9.78163185, 15.        , 21.32773109]), np.array([2879.68324683, 3227.54380282, 3612.83523974, 4030.93069323]), np.array([1453.20887391, 1924.88014667, 2383.65647311, 2842.30294301]))
        assert np.allclose(result[:1], expected_result[:1], rtol=1e-6, atol=1e-6)
        assert np.allclose(result[1:], expected_result[1:], rtol=1e-6, atol=1e-6)
    def test_Squirt_anisotropic(self):
        Sdry= [0.01,0,0.04,0.05,0.1]
        Sdry_hp= [0.015,0.1,0.045,0.15,0.15]
        result = Fluid.Squirt_anisotropic(Sdry, Sdry_hp)
        expected_result = np.array([0.07253425, 0.04520548, 0.04226027, 0.15547945, 0.16347603])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_White_Dutta_Ode(self,fluid_input):
        Kdry,Gdry,K0,Kfl,rho0,_,_,phi,_,a,_,freq= fluid_input
        rhofl1=1.05
        rhofl2=1.2
        Kfl1= 2.3
        Kfl2=2.4
        eta1= 0.001
        eta2= 0.002
        # kapa and a chosen so that |alpha*a| is about 0.2-20 over freq.
        # For very small alpha*a, Z1 and Z2 lose precision and the result is rounding noise.
        kapa= 1e-6
        a= 0.02
        sg= 0.5
        Vp, a_w, K_star = Fluid.White_Dutta_Ode(Kdry, Gdry, K0, phi, rho0, rhofl1,rhofl2, Kfl1, Kfl2,eta1,eta2,kapa,a,sg,freq )
        np.testing.assert_allclose(Vp, [2.838760946, 2.8387613368, 2.8388111118, 2.8388491456], rtol=1e-6)
        np.testing.assert_allclose(a_w, [1.9442463977e-08, 8.9735783941e-06, 1.3060960105e-03, 6.4526343269e-03], rtol=1e-6)
        np.testing.assert_allclose(K_star, [14.5846836597+3.4332042281e-06j, 14.5846890407+7.3549689390e-05j,
                                            14.5853743382+4.9691235595e-04j, 14.5858980091+1.1395302791e-04j], rtol=1e-6)
    def test_White_Dutta_Ode_low_frequency_limit(self,fluid_input):
        # Small alpha*a (about 1e-5 to 1e-3 here), where Z1 and Z2 used to lose precision.
        # The model must reduce to Gassmann with Wood's average of the two fluids,
        # and the attenuation must not be negative.
        Kdry,Gdry,K0,Kfl,rho0,_,_,phi,_,a,_,freq= fluid_input
        rhofl1, rhofl2, Kfl1, Kfl2 = 1.05, 1.2, 2.3, 2.4
        eta1, eta2, kapa, a, sg = 0.001, 0.002, 0.01, 0.0001, 0.5
        Vp, a_w, K_star = Fluid.White_Dutta_Ode(Kdry, Gdry, K0, phi, rho0, rhofl1,rhofl2, Kfl1, Kfl2,eta1,eta2,kapa,a,sg,freq )
        Kwood = 1/(sg/Kfl1 + (1-sg)/Kfl2)
        K_gassmann_wood = Kdry + (1-Kdry/K0)**2/(phi/Kwood + (1-phi)/K0 - Kdry/K0**2)
        np.testing.assert_allclose(K_star.real, K_gassmann_wood, rtol=1e-10)
        assert np.all(a_w >= 0)