#%%
#from import_path import rockphypy
from rockphypy import Anisotropy

import numpy as np


def _rotate_tensor(C, a):
    """Rotate a 6x6 Voigt stiffness matrix with the full tensor formula
    C'_ijkl = a_ip a_jq a_kr a_ls C_pqrs, independent of the Bond matrix."""
    voigt = {(0,0):0, (1,1):1, (2,2):2, (1,2):3, (2,1):3, (0,2):4, (2,0):4, (0,1):5, (1,0):5}
    T = np.zeros((3, 3, 3, 3))
    for (i, j), I in voigt.items():
        for (k, l), J in voigt.items():
            T[i, j, k, l] = C[I, J]
    T = np.einsum('ip,jq,kr,ls,pqrs->ijkl', a, a, a, a, T)
    C_rot = np.zeros((6, 6))
    for (i, j), I in voigt.items():
        for (k, l), J in voigt.items():
            C_rot[I, J] = T[i, j, k, l]
    return C_rot


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
        # a VTI medium rotated about its own symmetry axis (axis 3) must not change
        C =  np.array([[50, 30, 15,  0,  0,  0],
                    [30, 50, 15,  0,  0,  0],
                    [15, 15, 30,  0,  0,  0],
                    [ 0,  0,  0, 20,  0,  0],
                    [ 0,  0,  0,  0, 20,  0],
                    [ 0,  0,  0,  0,  0, 10]])
        result = Anisotropy.Bond_trans(C, 30, axis=3)
        np.testing.assert_allclose(result, C, atol=1e-12)
    def test_Bond_trans_tensor_rotation(self):
        # Bond transformation must give the same result as rotating the full
        # 4th-rank tensor, for a general anisotropic C, both axes and several angles
        rng = np.random.default_rng(0)
        A = rng.normal(size=(6, 6))
        C = 10*A@A.T + 50*np.eye(6)
        for theta in (17., 30., 90., 135.):
            t = np.deg2rad(theta)
            c, s = np.cos(t), np.sin(t)
            a3 = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])   # axis 3 fixed, same direction cosines as Bond_trans
            a1 = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])   # axis 1 fixed
            np.testing.assert_allclose(Anisotropy.Bond_trans(C, theta, axis=3), _rotate_tensor(C, a3), atol=1e-10)
            np.testing.assert_allclose(Anisotropy.Bond_trans(C, theta, axis=1), _rotate_tensor(C, a1), atol=1e-10)
