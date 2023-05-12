#%%
#from import_path import rockphypy
from rockphypy import AVO

import numpy as np



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
    def test_Aki_Richard(self):
        theta= np.linspace(0,40,10)
        # vp1,vp2,vs1,vs2= np.linspace(1500,2000,3), np.linspace(2500,3000,3), np.linspace(900,1000,3), np.linspace(1000,1200,3)
        # den1,den2= np.linspace(1.5,2.2,3),np.linspace(1.,2.2,3)
        vp1,vp2,vs1,vs2,den1,den2= 2000,3000,1000,1200,2200,2300
        result= AVO.Aki_Richard(theta, vp1,vp2,vs1,vs2,den1,den2)
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
        expected_result= (np.array([-0.33333333,  0.63570799,  0.2306779 ,  0.49624365,  0.62297421,
        -0.23233527,  0.58942721, -0.17987888,  0.558575  ,  0.46210946]), np.array([ 0.        ,  0.06837842,  0.3803706 ,  0.27702168, -0.05652195,
        -0.24909291, -0.1180928 , -0.29264611, -0.10438406,  0.21227566]))
        np.testing.assert_allclose(result, expected_result)
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
