#from import_path import rockphypy
from rockphypy import GM
import pytest
import numpy as np

class TestGM: 

    @pytest.fixture
    def GM_input(self):
        
        phic=0.4 # critical porosity
        sigma=20 # effective pressure 
        Cn=8.6 #coordination number
        f=0.5# reduced shear factor   
        K0, G0 = 36.6, 45 ## grain density, bulk and shear modulus 
        Kc, Gc = 36.6, 45 # cement density, bulk and shear modulus
        phib=0.3 # critical cement porosity
        return phic, sigma,Cn, f, K0, G0, Kc, Gc,phib
    def test_ThomasStieber(self):
        phi_sand=0.3
        phi_sh=0.2
        vsh=np.linspace(0,1,5)
        result = GM.ThomasStieber(phi_sand, phi_sh, vsh)
        expected_result =  (np.array([0.3 , 0.1 , 0.1 , 0.15, 0.2 ]),np.array([0.3  , 0.275, 0.25 , 0.225, 0.2  ]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    def test_contactcement(self,GM_input):

        phic, _,Cn, _, K0, G0, Kc, Gc,_ = GM_input
        phi = np.linspace(1e-7,phic,5) 
        result = GM.contactcement(K0, G0, Kc, Gc, phi, phic, Cn,  scheme=2)
        expected_result = (np.array([15.23060407, 13.38715662, 11.12410371,  8.05043534,  0.04975205]), np.array([20.63320692, 18.20079046, 15.18898382, 11.0578813 ,  0.12356857]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_hertzmindlin(self,GM_input):
        phic, sigma,Cn, f, K0, G0, Kc, Gc,phib = GM_input
        result = GM.hertzmindlin(K0, G0, phic, Cn, sigma, 1)
        expected_result = (1.9063198609211909, 2.802805417138183)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_softsand(self,GM_input):
        phic, sigma,Cn, f, K0, G0, Kc, Gc,phib = GM_input
        phi = np.linspace(1e-7,phic,5) 
        result = GM.softsand(K0, G0, phi, phic, Cn, sigma,f) 
        expected_result = (np.array([36.59992501, 10.84323519,  5.50294039,  3.19420155,  1.90631986]), np.array([44.99986823, 10.43728104,  5.22047177,  3.11303278,  1.97329867]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Walton(self,GM_input):
        phic, sigma,Cn, f, K0, G0, Kc, Gc,phib = GM_input
        phi = np.linspace(1e-7,phic,5) 
        result = GM.Walton(K0, G0, phic, Cn, sigma, f)
        expected_result = (1.9063198609211907, 1.9732986668454486)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_johnson(self,GM_input):
        K0,G0=37,44
        phi= 0.36
        Cn = 6
        epsilon= -0.03
        epsilon_axial = -0.05
        result = GM.johnson(K0, G0,Cn, phi, epsilon, epsilon_axial, path='together')
        expected_result = (np.array([[13.55173859,  2.99195167,  0.10545295,  0.        ,  0.        ,
          0.        ],
        [ 2.99195167, 13.55173859,  0.10545295,  0.        ,  0.        ,
          0.        ],
        [ 0.10545295,  0.10545295, 14.1114226 ,  0.        ,  0.        ,
          0.        ],
        [ 0.        ,  0.        ,  0.        ,  6.14910708,  0.        ,
          0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        ,  6.14910708,
          0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
          5.27989346]]), -0.7568273233538366, -0.21993424428637734)

        
        assert np.allclose(result[0], expected_result[0], rtol=1e-6, atol=1e-6)
        assert np.allclose(result[1:], expected_result[1:], rtol=1e-6, atol=1e-6)
    def test_stiffsand(self,GM_input):
        phic, sigma,Cn, f, K0, G0, Kc, Gc,phib = GM_input
        phi = np.linspace(1e-7,phic,5) 
        result = GM.stiffsand(K0, G0, phi, phic, Cn, sigma,f)
        expected_result = (np.array([36.59998647, 24.72898818, 15.45629949,  8.0129985 ,  1.90631986]), np.array([44.99997843, 27.76374785, 16.29187542,  8.1069418 ,  1.97329867]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_constantcement(self,GM_input):
        phic, sigma,Cn, f, K0, G0, Kc, Gc,phib = GM_input
        phi = np.linspace(1e-7,phic,5) 
        result = GM.constantcement(phib, K0, G0,Kc,Gc, phi, phic, Cn, scheme=2)
        expected_result = (np.array([36.59997856, 21.47765637, 13.23657067,  8.05043523,  4.48617491]), np.array([44.9999703 , 25.73592158, 16.48882383, 11.0578815 ,  7.48484988]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    # def test_MUHS(self,GM_input):
        phic, sigma,Cn, f, K0, G0, Kc, Gc,phib = GM_input
        phi = np.linspace(1e-7,phic,5) 
        result = GM.MUHS(K0, G0, Kc,Gc,phi, phib,phic, Cn, scheme=2)
        expected_result = (np.array([36.59998649, 24.7483484 , 15.48701406,  8.05043461,  1.94767854]), np.array([44.99998129, 29.63592467, 18.93672859, 11.05788087,  5.01406624]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Digby(self):
        a_R=0.2
        G0= 44
        K0=37
        D0= 2.65
        phi= 0.36
        D= (1-phi)* D0
        Cn=9
        sigma=20# Mpa
        result = GM.Digby(K0, G0, phi, Cn, sigma, a_R )
        expected_result = (np.array([5.81155426]), np.array([8.51375906]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_pcm(self,GM_input):
        phic, sigma,Cn, f, K0, G0, Kc, Gc,phib = GM_input
        phi=np.linspace(1e-7,phic,5)
        v_cem=0.1
        v_ci=0.111
        scheme=2
        f_=0 #reduce shear factor
        result = GM.pcm(f,sigma, K0,G0,phi, phic,v_cem,v_ci, Kc,Gc,Cn=Cn, mode='stiff',scheme=scheme,f_=f_)
        expected_result = (np.array([36.59996762, 18.18928411, 10.87569332,  6.94950669,  4.49985482]), np.array([44.9999437 , 18.71837729, 10.71568312,  6.84274853,  4.55850738]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)