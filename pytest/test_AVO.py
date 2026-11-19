#%%
#from import_path import rockphypy
from rockphypy import AVO

import numpy as np
import pytest





def _zoeppritz_boundary_solution(vp1, vs1, rho1, vp2, vs2, rho2, theta):
    """Solve displacement/traction continuity for a unit incident P wave.

    Independent of the closed-form a/b/c/d/E/F/G/H implementation.
    Matrix and sign convention: rockphypy's AVO theory documentation,
    https://rockphypy.readthedocs.io/en/latest/getting_started/11_AVO.html
    (M and the first column of N, attributed there to Aki & Richards, 1980).
    Only real, subcritical solid/solid interfaces are used in these tests.
    """
    theta1 = np.deg2rad(theta)
    p = np.sin(theta1) / vp1
    assert max(abs(p * vp2), abs(p * vs1), abs(p * vs2)) < 1
    theta2, phi1, phi2 = np.arcsin([p * vp2, p * vs1, p * vs2])
    s1, c1 = np.sin(theta1), np.cos(theta1)
    s2, c2 = np.sin(theta2), np.cos(theta2)
    q1, d1 = np.sin(phi1), np.cos(phi1)
    q2, d2 = np.sin(phi2), np.cos(phi2)
    matrix = np.array([
        [-s1, -d1, s2, d2],
        [c1, -q1, c2, -q2],
        [2*rho1*vs1*q1*c1, rho1*vs1*(1-2*q1*q1),
         2*rho2*vs2*q2*c2, rho2*vs2*(1-2*q2*q2)],
        [-rho1*vp1*(1-2*q1*q1), 2*rho1*vs1*q1*d1,
         rho2*vp2*(1-2*q2*q2), -2*rho2*vs2*q2*d2],
    ])
    rhs = np.array([s1, c1, 2*rho1*vs1*q1*c1,
                    rho1*vp1*(1-2*q1*q1)])
    # Balance stress and displacement rows before solving.
    scale = np.max(np.abs(matrix), axis=1)
    return np.linalg.solve(matrix / scale[:, None], rhs / scale)


class TestAVO: 
    def test_AVO_HTI(self):
        C = np.array([[100, 30, 15,  0,  0,  0],
              [30, 100, 15,  0,  0,  0],
              [15, 15, 30,  0,  0,  0],
              [0,  0,  0, 20,  0,  0],
              [0,  0,  0,  0, 20,  0],
              [0,  0,  0,  0,  0, 10]])
        D = np.array([[97.,  7.,  7.,  0.,  0.,  0.],
                    [7., 97.,  7.,  0.,  0.,  0.],
                    [7.,  7., 97.,  0.,  0.,  0.],
                    [0.,  0.,  0., 45.,  0.,  0.],
                    [0.,  0.,  0.,  0., 45.,  0.],
                    [0.,  0.,  0.,  0.,  0., 45.]])
        theta= np.linspace(0,20,4)
        azimuth= np.linspace(0,20,4)

        result = AVO.AVO_HTI(2.2,2.4,D,C,theta,azimuth)
        expected_result = np.array(
        [[-0.26513677, -0.23127679, -0.13089322,  0.03261608],
        [-0.26513677, -0.23164028, -0.13233472,  0.02941735],
        [-0.26513677, -0.2327113 , -0.13658373,  0.01998204],
        [-0.26513677, -0.23443248, -0.14341743,  0.00478649]])
        np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Aki_Richards(self):
        theta= np.linspace(0,40,10)
        # vp1,vp2,vs1,vs2= np.linspace(1500,2000,3), np.linspace(2500,3000,3), np.linspace(900,1000,3), np.linspace(1000,1200,3)
        # den1,den2= np.linspace(1.5,2.2,3),np.linspace(1.,2.2,3)
        vp1,vp2,vs1,vs2,den1,den2= 2000,3000,1000,1200,2200,2300
        result= AVO.Aki_Richards(theta, vp1,vp2,vs1,vs2,den1,den2)
        expected_result= (np.array([0.22222222, 0.22248164, 0.22334161, 0.22505294, 0.22805364,
        0.2330032 , 0.24084072, 0.25288048, 0.27096972, 0.29775439]),np.array([-0.        , -0.01942277, -0.03799843, -0.05490899, -0.069394  ,
        -0.08077786, -0.08849507, -0.09211298, -0.09135103, -0.08609571]), 0.22222222222222224,
        0.04199111111111112) 
        # 
        np.allclose(result[:2], expected_result[:2], rtol=1e-6, atol=1e-6)
        np.allclose(result[2:], expected_result[2:], rtol=1e-6, atol=1e-6)
    def test_zoeppritz(self):
        theta= np.linspace(0,40,10)
        vp1, vs1, rho1, vp2, vs2, rho2, = 4400,3200,2500, 2500,1300,2200
        result= AVO.zoeppritz(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        expected_result = np.array([
            _zoeppritz_boundary_solution(vp1, vs1, rho1, vp2, vs2, rho2, angle)[:2]
            for angle in theta
        ]).T
        np.testing.assert_allclose(result, expected_result, rtol=1e-12, atol=1e-12)

    @pytest.mark.parametrize("properties", [
        (2400., 1100., 2150., 3100., 1700., 2450.),
        (3600., 2000., 2500., 2700., 1400., 2200.),
        (3000., 1600., 2300., 3000., 1200., 2300.),
    ])
    @pytest.mark.parametrize("angle", [0., 7., 20., 35.])
    def test_zoeppritz_boundary_conditions(self, properties, angle):
        expected = _zoeppritz_boundary_solution(*properties, angle)
        actual = AVO.zoeppritz(*properties, angle)
        np.testing.assert_allclose(actual, expected[:2], rtol=1e-12, atol=1e-12)

    def test_zoeppritz_physical_limits(self):
        properties = (2400., 1100., 2150., 3100., 1700., 2450.)
        rpp, rps = AVO.zoeppritz(*properties, 0.)
        z1, z2 = properties[0] * properties[2], properties[3] * properties[5]
        np.testing.assert_allclose([rpp, rps], [(z2-z1)/(z2+z1), 0.], atol=1e-12)
        theta = np.array([0., 7., 20., 35.])
        np.testing.assert_allclose(
            AVO.zoeppritz(3000., 1600., 2300., 3000., 1600., 2300., theta),
            np.zeros((2, theta.size)), atol=1e-12)

    def test_zoeppritz_angle_symmetry(self):
        properties = (2400., 1100., 2150., 3100., 1700., 2450.)
        theta = np.array([7., 20., 35.])
        positive = AVO.zoeppritz(*properties, theta)
        negative = AVO.zoeppritz(*properties, -theta)
        np.testing.assert_allclose(positive[0], negative[0], atol=1e-12)
        np.testing.assert_allclose(positive[1], -negative[1], atol=1e-12)

    def test_AVO_abe(self):
        vp1,vp2,vs1,vs2= np.linspace(1500,2000,3), np.linspace(2500,3000,3), np.linspace(900,1000,3), np.linspace(1000,1200,3)
        den1,den2= np.linspace(1.5,2.2,3),np.linspace(1.,2.2,3)
        result= AVO.AVO_abe(vp1,vs1,den1,vp2,vs2,den2)
        expected_result = (np.array([0.05, 0.14975845, 0.2]),
                           np.array([0.31447918, 0.14255428, 0.04734781]),
                           np.array([0.3355, 0.16089461, 0.0592]),
                           np.array([0.27984688,  0.01213751, -0.144512]),
                           np.array([0.36631579,  0.00614048, -0.18181818]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    def test_EI_ref(self):
        theta= np.linspace(0,40,4)
        Vp, Vs= np.linspace(2500,3000,3),np.linspace(1500,2000,3)
        SP=0.5
        rho= np.linspace(2.000,2.200,3)
        result= AVO.EI_ref(Vp,Vs,rho,theta,SP,norm=True)
        expected_result = (np.array([[5000.        , 5068.65592815, 5245.25596794, 5418.791153  ],
        [5775.        , 5775.        , 5775.        , 5775.        ],
        [6600.        , 6522.64313991, 6333.50579279, 6164.26231533]]), np.array([[5775.        , 6314.14415614, 6703.71276525, 6769.41465685],
        [5775.        , 5775.        , 5775.        , 5775.        ],
        [5775.        , 5335.09168518, 5057.60056572, 5010.44840975]]), np.array([[3675.        , 3878.87409752, 4027.7227177 , 4030.1759451 ],
        [3675.        , 3675.        , 3675.        , 3675.        ],
        [3675.        , 3486.13535395, 3372.38883288, 3446.06997545]]), np.array([[4501.875     , 4393.62795852, 4108.14932587, 3739.56534509],
        [3675.        , 3675.        , 3675.        , 3675.        ],
        [3069.46022727, 3151.62398647, 3389.42657025, 3750.17554161]]), np.array([[4501.875     , 4488.3492232 , 4447.94597164, 4382.23399945],
        [3675.        , 3675.        , 3675.        , 3675.        ],
        [3069.46022727, 3079.39898699, 3110.02516231, 3163.05498124]]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
        
#%%    
