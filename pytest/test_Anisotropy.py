#%%
#from import_path import rockphypy
from rockphypy import Anisotropy

import numpy as np

#%%
class TestAnisotropy:
    def test_Thomsen(self):
        C11, C33, C13, C44, C66, den, theta = 20, 20, 10, 10, 20, 2, 30.0
        expected_result = (3606.972956129558, 1607.1738588279738, 2515.5764746872637, 0.0, 0.5, 0.75)
        result = Anisotropy.Thomsen(C11, C33, C13, C44, C66, den, theta)
        np.testing.assert_allclose(result, expected_result)

    def test_Thomsen_Tsvankin(self):
        C11,C22,C33,C12,C13,C23,C44,C55,C66 = 20, 5,20, 10, 10, 20, 25, 30,15 
        expected_result = (-0.375, -10.0, -0.25, 0.0, -3.75, -0.2, 3.0)
        result = Anisotropy.Thomsen_Tsvankin(C11,C22,C33,C12,C13,C23,C44,C55,C66)
        np.testing.assert_allclose(result, expected_result)
    def test_Backus(self):
        V=[0.1,0.2,0.7]
        lamda=[15,20,10]
        G=[20,10,15]
        result = Anisotropy.Backus(V, lamda, G)
        expected_result= (41.098130841121495,
                        41.12149532710281,
                        12.429906542056075,
                        13.953488372093023,
                        14.5)
        np.testing.assert_allclose(result, expected_result)
    def test_Backus_log(self):
        
        Vp=[3.000,3.500,4.000]
        Vs=[1.550,1.320,1.400]
        Den= [2.2,2.3,2.35]
        Depth=[1000,1700,2300,4000]
        result = Anisotropy.Backus_log(Vp,Vs,Den,Depth)
        expected_result= (29.133954660361283,
                            29.45167576239999,
                            20.018242129339992,
                            4.6065956954832,
                            4.6448540000000005,
                            2.305)
        np.testing.assert_allclose(result, expected_result)
    def test_vel_azi_HTI(self):
        
        C=np.array([[30, 15, 15,  0,  0,  0],
                    [15, 30, 10,  0,  0,  0],
                    [15, 10, 30,  0,  0,  0],
                    [ 0,  0,  0, 10,  0,  0],
                    [ 0,  0,  0,  0, 25,  0],
                    [ 0,  0,  0,  0,  0, 25]])
        Den= 2.3
        azimuth= np.linspace(0,90,10)
        result = Anisotropy.vel_azi_HTI(C,Den,azimuth)
        expected_result= (np.array([3.61157559, 3.88601535, 
                4.19614812, 4.41664613, 4.53013766,
                4.53013766, 4.41664613, 4.19614812, 3.88601535, 3.61157559]), np.array([3.29690237, 2.9684892 , 2.51105246, 2.0991142 , 1.84143863,
                1.84143863, 2.0991142 , 2.51105246, 2.9684892 , 3.29690237]), np.array([3.29690237, 3.26694211, 3.1790984 , 3.03959379, 2.85918515,
                2.65376178, 2.44504823, 2.26069119, 2.13177874, 2.08514414]))
        np.testing.assert_allclose(result, expected_result)
    def test_vel_azi_VTI(self):
        C =  np.array([[50, 30, 15,  0,  0,  0],
                    [30, 50, 15,  0,  0,  0],
                    [15, 15, 30,  0,  0,  0],
                    [ 0,  0,  0, 20,  0,  0],
                    [ 0,  0,  0,  0, 20,  0],
                    [ 0,  0,  0,  0,  0, 10]])
        Den= 2.3
        azimuth= np.linspace(0,90,10)
        result = Anisotropy.vel_azi_VTI(C,Den,azimuth)
        expected_result= (np.array([3.61157559, 3.77222919, 
                4.05016857, 4.30648876, 4.50469172,
                4.63315913, 4.69307069, 4.69810515, 4.67602567, 4.66252404]), np.array([2.94883912, 2.9265251 , 2.8613027 , 2.75838642, 2.62663952,
                2.47875838, 2.33126202, 2.20372982, 2.11634805, 2.08514414]), np.array([2.94883912, 2.78776315, 2.5204093 , 2.31672139, 2.2449314 ,
                2.31857393, 2.49718983, 2.71023847, 2.88224921, 2.94883912]))    
        np.testing.assert_allclose(result, expected_result)   
    def test_Bond_trans(self):
        C =  np.array([[50, 30, 15,  0,  0,  0],
                    [30, 50, 15,  0,  0,  0],
                    [15, 15, 30,  0,  0,  0],
                    [ 0,  0,  0, 20,  0,  0],
                    [ 0,  0,  0,  0, 20,  0],
                    [ 0,  0,  0,  0,  0, 10]])
        theta=30
        result = Anisotropy.Bond_trans(C, theta, axis=3)
        expected_result=np.array([[ 4.43750000e+01,  3.56250000e+01,  1.50000000e+01,
          0.00000000e+00,  0.00000000e+00, -2.16506351e+00],
        [ 3.56250000e+01,  4.43750000e+01,  1.50000000e+01,
          0.00000000e+00,  0.00000000e+00,  2.16506351e+00],
        [ 1.50000000e+01,  1.50000000e+01,  3.00000000e+01,
          0.00000000e+00,  0.00000000e+00,  3.88578059e-16],
        [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
          2.00000000e+01,  1.52242806e-15,  0.00000000e+00],
        [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
          2.06956853e-15,  2.00000000e+01,  0.00000000e+00],
        [-2.16506351e+00,  2.16506351e+00,  3.88578059e-16,
          0.00000000e+00,  0.00000000e+00,  1.00000000e+01]])
        np.testing.assert_allclose(result, expected_result)  
