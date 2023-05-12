
#%%
#from import_path import rockphypy
from rockphypy import utils
import pytest
import numpy as np


#%%
class Testutils:

    @pytest.fixture
    def input_mod(self):
        K,G= 37,44
        rho=2.65
        return K,G, rho 

    @pytest.fixture
    def input_C(self):
        C11,C22,C33,C12,C13,C23,C44,C55,C66 = 40, 5,20, 10, 10, 20, 5, 30,15 
        return C11,C22,C33,C12,C13,C23,C44,C55,C66
    def test_V(self,input_mod):
        K,G, rho = input_mod
        result = utils.V(K, G, rho)
        expected_result = (6008.379892351815, 4074.7728261714988)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_poi(self,input_mod):
        K,G, _ = input_mod
        result = utils.poi(K, G)
        expected_result = 0.07419354838709677
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    def test_lame(self,input_mod):
        K,G, _ = input_mod
        result = utils.lame(K, G)
        expected_result = 7.666666666666668
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_M_from_V(self,input_mod):
        _,_, rho = input_mod
        result = utils.M_from_V(rho, 6000,4000)
        expected_result = (38.866666666666674, 42.4)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_write_HTI_matrix(self, input_C):
        C11,C22,C33,C12,C13,C23,C44,C55,C66 =  input_C
        result = utils.write_HTI_matrix(C11,C33,C13,C44,C55)
        expected_result = np.array([[40, 10, 10,  0,  0,  0],
                                    [10, 20, 10,  0,  0,  0],
                                    [10, 10, 20,  0,  0,  0],
                                    [ 0,  0,  0,  5,  0,  0],
                                    [ 0,  0,  0,  0, 30,  0],
                                    [ 0,  0,  0,  0,  0, 30]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    def test_write_VTI_compliance(self,input_C):
        C = np.array([[40, 10, 10,  0,  0,  0],
                    [10, 40, 10,  0,  0,  0],
                    [10, 10, 20,  0,  0,  0],
                    [ 0,  0,  0,  5,  0,  0],
                    [ 0,  0,  0,  0,  5,  0],
                    [ 0,  0,  0,  0,  0, 15]])
        S = np.linalg.inv(C)
        result = utils.write_VTI_compliance(S[0,0],S[0,1],S[0,2],S[2,2],S[3,3])
        expected_result = np.array([[ 0.02916667, -0.00416667, -0.0125    ,  0.        ,  0.        ,
         0.        ],
       [-0.00416667,  0.02916667, -0.0125    ,  0.        ,  0.        ,
         0.        ],
       [-0.0125    , -0.0125    ,  0.0625    ,  0.        ,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.2       ,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.2       ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         0.06666667]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    def test_write_VTI_matrix(self,input_C):
        C11,C22,C33,C12,C13,C23,C44,C55,C66 =  input_C
        result = utils.write_VTI_matrix(C11,C33,C13,C44,C66)
        expected_result = np.array([[40, 10, 10,  0,  0,  0],
                                    [10, 40, 10,  0,  0,  0],
                                    [10, 10, 20,  0,  0,  0],
                                    [ 0,  0,  0,  5,  0,  0],
                                    [ 0,  0,  0,  0,  5,  0],
                                    [ 0,  0,  0,  0,  0, 15]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    
    def test_write_matrix(self,input_C):
        C11,C22,C33,C12,C13,C23,C44,C55,C66 =  input_C
        result = utils.write_matrix(C11,C22,C33,C12,C13,C23,C44,C55,C66)
        expected_result = np.array([[40, 10, 10,  0,  0,  0],
                                    [10,  5, 20,  0,  0,  0],
                                    [10, 20, 20,  0,  0,  0],
                                    [ 0,  0,  0,  5,  0,  0],
                                    [ 0,  0,  0,  0, 30,  0],
                                    [ 0,  0,  0,  0,  0, 15]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_write_iso(self,input_mod):
        K,G, _ = input_mod
        result = utils.write_iso(K,G)
        expected_result = np.array([[95.66666667,  7.66666667,  7.66666667,  0.        ,  0.        ,
         0.        ],
       [ 7.66666667, 95.66666667,  7.66666667,  0.        ,  0.        ,
         0.        ],
       [ 7.66666667,  7.66666667, 95.66666667,  0.        ,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        , 44.        ,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        , 44.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        44.        ]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    def test_crack_por(self):
        result = utils.crack_por(0.01, 0.1)
        expected_result = 0.004188790204786391
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    def test_v_to_c_VTI(self):
        result = utils.v_to_c_VTI(5.500,2.800,4.000,3.000,2.800,2.5)
        expected_result = np.array([[40.        ,  0.8       , 14.55337367,  0.        ,  0.        ,
          0.        ],
        [ 0.8       , 40.        , 14.55337367,  0.        ,  0.        ,
          0.        ],
        [14.55337367, 14.55337367, 75.625     ,  0.        ,  0.        ,
          0.        ],
        [ 0.        ,  0.        ,  0.        , 22.5       ,  0.        ,
          0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        , 22.5       ,
          0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         19.6       ]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    


# %%
