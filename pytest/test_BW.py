#from import_path import rockphypy
from rockphypy import BW

import numpy as np
import pytest

class TestBW:

    @pytest.fixture
    def input_data(self):
        P= np.linspace(20,40,4)
        T= np.linspace(40,100,4)
        G= 1.2
        den = 0.88
        return P, T, G,den
    def test_rho_K_co2(self, input_data):
        P, T, G,_  = input_data
        result = BW.rho_K_co2(P,T,G)
        expected_result= (np.array([0.61592007, 0.60874166, 0.60099684, 0.59403165]), np.array([0.11058432, 0.15243223, 0.18276998, 0.20478121]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    def test_rho_K_gas(self,input_data):
        P, T, G,_  = input_data
        result = BW.rho_K_gas(P,T,G)
        expected_result= (np.array([0.41957204, 0.41635711, 0.41348241, 0.41128399]), np.array([0.16693798, 0.2037523 , 0.22983503, 0.24951797]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_rho_K_oil(self,input_data): 
        P, T,_ ,den  = input_data
        result = BW.rho_K_oil(P,T,den)
        expected_result= (np.array([0.87622229, 0.86380022, 0.85103238, 0.83806654]),np.array([1.84963494, 1.73858842, 1.63828494, 1.54822048]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_rho_K_go(self,input_data):
        P, T,G ,den  = input_data
        Rg= 85
        result = BW.rho_K_go(P,T,den,G,Rg)
        expected_result= (np.array([0.78527236, 0.77059486, 0.75630133, 0.74238677]),np.array([1.05535966, 0.97603187, 0.91739857, 0.8784189 ]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_rho_K_water(self,input_data): 
        P, T,_,_ = input_data
        result = BW.rho_K_water(T,P)
        expected_result= (np.array([1.00014216, 0.99467699, 0.98698636, 0.9773372 ]),np.array([2.4395262 , 2.5457876 , 2.59353244, 2.59199095]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_v_water(self,input_data):
        P, T,_,_ = input_data
        result = BW.v_water(T,P)
        expected_result= (np.array([1561.78726157, 1599.81604506, 1621.0270707 , 1628.52536   ]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    def test_rho_K_brine(self,input_data):
        S= 240000/1000000
        P, T,_,_ = input_data
        result = BW.rho_K_brine(T,P,S)
        expected_result= (np.array([1.17846792, 1.17025331, 1.16027916, 1.1488124 ]), np.array([3.84408197, 3.87376277, 3.86013756, 3.8135643 ]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_v_brine(self,input_data):
        S= 240000/1000000
        P, T,_,_ = input_data
        result = BW.v_brine(T,P,S)
        expected_result= (np.array([1806.0818885 , 1819.39321897, 1823.98030992, 1821.96898287]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)


    