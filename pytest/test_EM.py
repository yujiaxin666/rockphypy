#from import_path import rockphypy
from rockphypy import EM

import numpy as np
from scipy.integrate import solve_ivp
from rockphypy.utils import utils


class TestEM:
    def test_VRH(self): 
        volumes=[0.1,0.2,0.7]
        M= [10,20,30]
        result=EM.VRH(volumes, M)
        expected_result = (26.0, 23.076923076923077, 24.53846153846154)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_cripor(self):
        K0, G0= 47,44
        phic=0.4
        phi= np.linspace(1e-6,phic,10 )
        result= EM.cripor(K0, G0, phi, phic)
        expected_result = (np.array([46.9998825 , 41.77767333, 36.55546417, 31.333255  , 26.11104583, 20.88883667, 15.6666275 , 10.44441833,  5.22220917,  0.        ]),
        np.array([43.99989   , 39.11101333, 34.22213667, 29.33326   , 24.44438333, 19.55550667, 14.66663   ,  9.77775333,  4.88887667,  0.        ]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_cripor_reuss(self):
        phic=0.4
        M0= 10
        Mf= 2.2
        result = EM.cripor_reuss(M0, Mf, phic)
        expected_result = np.array([4.13533835])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)


    def test_HS(self): 
        f= 0.5
        K1=np.linspace(40,50,2)
        K2=np.linspace(35,37,2)
        G1=np.linspace(20,30,2)
        G2=np.linspace(15,20,2)
        result = EM.HS(f, K1, K2,G1, G2, bound='upper')
        expected_result = (np.array([37.4025974 , 42.99401198]), np.array([17.34042553, 24.55645161]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Eshelby_Cheng(self):
        K, G= 37,45
        phi=0.01
        alpha= 0.1
        Kf= 2.1
        result = EM.Eshelby_Cheng(K, G, phi, alpha, Kf, mat=False)
        expected_result = np.array([96.00459425, 86.77207515,  6.82580934, 41.68262752, 44.48956533])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_hudson(self):
        K, G= 37,45
        Ki, Gi = 2.1, 0
        alpha= 0.1
        crd = 0.01
        result = EM.hudson(K, G, Ki, Gi, alpha, crd, order=1, axis=3)
        expected_result = np.array(
            [[96.97879354,  6.97879354,  6.70613902, 0., 0., 0.],
             [6.97879354, 96.97879354,  6.70613902,  0.,  0., 0.],
             [6.70613902,  6.70613902, 92.92792646,  0.,  0., 0.],
             [0.,  0.,  0., 43.84179104,  0., 0.],
             [0.,  0.,  0.,  0., 43.84179104, 0.],
             [0.,  0.,  0.,  0.,  0., 45.]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_hudson_rand(self):
        K, G= 37,45
        Ki, Gi = 2.1, 0
        alpha= 0.1
        crd = 0.01
        result = EM.hudson_rand(K, G, Ki, Gi, alpha, crd)
        expected_result = (36.407517410762416, 44.30301254779482)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_hudson_rand_orientation_average(self):
        # First-order corrections are linear in crack density, so randomly oriented
        # cracks must give the aligned correction (EM.hudson) averaged over all orientations.
        K, G, alpha, crd = 37., 44., 0.01, 0.01
        voigt = {(0,0):0, (1,1):1, (2,2):2, (1,2):3, (2,1):3, (0,2):4, (2,0):4, (0,1):5, (1,0):5}
        for Ki, Gi in ((0., 0.), (2.25, 0.), (5., 3.)):
            C1 = EM.hudson(K, G, Ki, Gi, alpha, crd, order=1, axis=3) - utils.write_iso(K, G)
            T = np.zeros((3, 3, 3, 3))
            for (i, j), I in voigt.items():
                for (k, l), J in voigt.items():
                    T[i, j, k, l] = C1[I, J]
            # orientation average of a 4th-rank tensor = its isotropic part
            iijj, ijij = np.einsum('iijj', T), np.einsum('ijij', T)
            dG = (3*ijij - iijj)/30
            dK = (2*iijj - ijij)/15 + 2*dG/3
            K_eff, G_eff = EM.hudson_rand(K, G, Ki, Gi, alpha, crd)
            np.testing.assert_allclose([K_eff - K, G_eff - G], [dK, dG], rtol=1e-10)
    def test_hudson_ortho(self):
        K, G= 37,45
        Ki, Gi = 2.1, 0
        alpha= 0.1
        crd = 0.01
        alpha= [0.1,0.01,1]
        crd=[0.001,0.01,0.1]
        result= EM.hudson_ortho(K, G, Ki, Gi, alpha, crd)
        expected_result = np.array([[96.32208655,  6.60756844,  3.31174022,  0.        ,  0.        ,
         0.        ],
       [ 6.60756844, 95.35527413,  3.24666631,  0.        ,  0.        ,
         0.        ],
       [ 3.31174022,  3.24666631, 46.38868342,  0.        ,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        , 32.25970149,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        , 33.30208955,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        43.72597015]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_hudson_cone(self):
        K, G= 37,45
        Ki, Gi = 2.1, 0
        crd = 0.01
        alpha= 0.1
        #alpha= [0.1,0.01,1]
        crd=0.01#[0.001,0.01,0.1]
        theta = 30
        result= EM.hudson_cone(K, G, Ki, Gi, alpha, crd, theta)
        expected_result = np.array([[96.357946  ,  6.91943677,  6.84590313,  0.        ,  0.        ,
         0.        ],
       [ 6.91943677, 96.357946  ,  6.84590313,  0.        ,  0.        ,
         0.        ],
       [ 6.84590313,  6.84590313, 93.72927864,  0.        ,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        , 44.09224946,  0.        ,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        , 44.09224946,
         0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        44.6370931 ]])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_hudson_filled_cracks(self):
        # cracks filled with a weak solid (Gi > 0), so the M term matters
        K, G, Ki, Gi = 37., 44., 5., 3.
        alpha, crd = 0.01, 0.05
        lam = K - 2*G/3
        # M and U1 from Hudson (1981), Rock Physics Handbook section 4.10
        M = 4*Gi*(lam + 2*G)/(np.pi*alpha*G*(3*lam + 4*G))
        U1 = 16*(lam + 2*G)/(3*(3*lam + 4*G)*(1 + M))
        C = EM.hudson(K, G, Ki, Gi, alpha, crd, order=1, axis=3)
        np.testing.assert_allclose(C[3,3], G - G*crd*U1, rtol=1e-12)
        # one aligned crack set: hudson_cone (theta=0) and hudson_ortho must agree with hudson
        np.testing.assert_allclose(EM.hudson_cone(K, G, Ki, Gi, alpha, crd, 0), C, atol=1e-12)
        np.testing.assert_allclose(EM.hudson_ortho(K, G, Ki, Gi, np.array([alpha]*3), np.array([0, 0, crd])), C, atol=1e-12)
    def test_Berryman_sc(self):
        K= [47,35]
        G= [30,20]
        X= [0.3,0.7]
        Alpha = [0.1,1.2]
        result= EM.Berryman_sc(K,G,X,Alpha)
        expected_result = (38.20950220104562, 22.605875876194837)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_PQ_prolate_limits(self):
        Km, Gm, Ki, Gi = 37., 44., 2.25, 0.
        # P and Q must be continuous across alpha = 1 (sphere)
        P1, Q1 = EM.PQ(Km, Gm, Ki, Gi, 1)
        np.testing.assert_allclose(EM.PQ(Km, Gm, Ki, Gi, 0.999), (P1, Q1), rtol=1e-5)
        np.testing.assert_allclose(EM.PQ(Km, Gm, Ki, Gi, 1.001), (P1, Q1), rtol=1e-5)
        # very long prolate spheroids -> needle formulas (Mavko et al., Rock Physics Handbook)
        gam = Gm*(3*Km + Gm)/(3*Km + 7*Gm)
        P_needle = (Km + Gm + Gi/3)/(Ki + Gm + Gi/3)
        Q_needle = (4*Gm/(Gm + Gi) + 2*(Gm + gam)/(Gi + gam) + (Ki + 4*Gm/3)/(Ki + Gm + Gi/3))/5
        np.testing.assert_allclose(EM.PQ(Km, Gm, Ki, Gi, 1e4), (P_needle, Q_needle), rtol=1e-6)
    def test_PQ_vectorize(self):
        Km, Gm = 37,44
        Ki= np.linspace(2,3,4)
        Gi= np.zeros(4)
        alpha= np.array([0.005,0.001,0.01,0.1])
        result= EM.PQ_vectorize(Km,Gm, Ki,Gi, alpha)
        expected_result = (np.array([15.72732627, 15.39575083, 11.02069818,  3.90846602]), np.array([ 55.03822226, 250.79531002,  29.06865426,   4.82746831]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_PQ_vectorize_does_not_change_input(self):
        Km, Gm = 37., 44.
        Ki, Gi = np.zeros(2), np.zeros(2)
        for alpha in (np.array([1.0, 0.1]), np.array([1, 2])):
            alpha_in = alpha.copy()
            P, Q = EM.PQ_vectorize(Km, Gm, Ki, Gi, alpha)
            np.testing.assert_array_equal(alpha, alpha_in)                      # input unchanged
            np.testing.assert_allclose((P[0], Q[0]), EM.PQ(Km, Gm, 0., 0., 1))  # exact sphere values
    def test_Swiss_cheese(self):
        Ks=np.array([30,40,50])
        Gs=np.array([10,20,30])
        phi= np.linspace(0.1,0.3,5)
        result= EM.Swiss_cheese(Ks,Gs,phi)
        expected_result = (np.array([[22.64150943, 20.16806723, 18.18181818, 16.55172414, 15.18987342],
        [32.        , 29.09090909, 26.66666667, 24.61538462, 22.85714286],
        [40.81632653, 37.38317757, 34.48275862, 32.        , 29.85074627]]), np.array([[ 8.43373494,  7.82122905,  7.29166667,  6.82926829,  6.42201835],
        [16.77419355, 15.52238806, 14.44444444, 13.50649351, 12.68292683],
        [25.09090909, 23.19327731, 21.5625    , 20.1459854 , 18.90410959]]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_SC(self):
        phi= np.linspace(0.1,0.3,3)
        Ks, Gs = 37,44
        iter_n = 10
        result = EM.SC(phi,Ks,Gs,iter_n)
        expected_result =(np.array([30.84267143, 24.35677881, 17.39208156]), np.array([34.82985709, 25.77929529, 16.94540536]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Dilute_crack(self):
        Ks, Gs = 37,44
        cd= np.array([0.001,0.01,0.1])
        result = EM.Dilute_crack(Ks,Gs,cd)
        expected_result =(np.array([36.92334526, 36.24748381, 30.63915556]),np.array([43.92603192, 43.27134365, 37.65858256]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_OConnell_Budiansky(self):
        Ks, Gs = 37,44
        cd= np.array([0.001,0.01,0.1])
        result= EM.OConnell_Budiansky(Ks,Gs,cd)
        expected_result =(np.array([36.9232084 , 36.2340839 , 29.53604499]), np.array([43.9258999 , 43.25832733, 36.51631015]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_OConnell_Budiansky_fl(self):
        Ks, Gs = 37,44
        Kfl= 2.1
        crd= 0.1
        alpha= 0.1
        result= EM.OConnell_Budiansky_fl(Ks,Gs,Kfl,crd, alpha)
        expected_result = (31.248892685984856, 37.30710074731592)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)

    def test_PQ(self):
        Km, Gm = 37,44
        Ki, Gi = 2.1,0
        alpha = 0.1
        result = EM.PQ(Km,Gm, Ki,Gi, alpha)
        expected_result = (4.234473045076649, 4.924519092069241)
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Berryman_DEM(self):
        Gi =0
        Ki =0 # dry pore
        alpha = 1
        Km, Gm = 37,44
        phi =  0.1#np.arange(0,phic,0.001)
        result = EM.Berryman_DEM(Km,Gm, Ki, Gi, alpha,phi)
        # values checked against an independent ODE solve, see test_Berryman_DEM_sphere_reference
        expected_result = (np.array([37.          , 36.3980156421, 35.7987912512, 35.2023676466,
        34.6087847931, 34.0180832062, 33.4303041747, 32.8454896867,
        32.2636823319, 31.6849253783, 31.1092627245]), np.array([44.          , 43.0835020317, 42.1774950805, 41.2819511623,
        40.3968433271, 39.5221446883, 38.6578282172, 37.8038668939,
        36.9602337064, 36.1269016502, 35.3038437337]), np.array([0.  , 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 ]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Berryman_DEM_sphere_reference(self):
        # Independent check: integrate the DEM equations for dry spherical pores
        # with solve_ivp, using the closed-form P and Q for spheres
        # (not EM.PQ or EM.DEM).
        Km, Gm, Ki, Gi, phi = 37., 44., 0., 0., 0.1
        def rhs(t, y):
            K, G = y
            P = (K + 4*G/3) / (Ki + 4*G/3)
            zeta = G/6 * (9*K + 8*G) / (K + 2*G)
            Q = (G + zeta) / (Gi + zeta)
            return [(Ki - K)*P/(1 - t), (Gi - G)*Q/(1 - t)]
        K, G, t = EM.Berryman_DEM(Km, Gm, Ki, Gi, 1, phi)
        ref = solve_ivp(rhs, (0, t[-1]), [Km, Gm], t_eval=t, rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(K, ref.y[0], rtol=1e-6)
        np.testing.assert_allclose(G, ref.y[1], rtol=1e-6)
    def test_SC_dilute(self):
        ro=13/6
        Km, Gm=3*ro, 3 #
        Ki, Gi=Km*50, Gm*50  #  65, 30
        f= np.linspace(0,0.4,5)
        #  Dilute distribution of inclusion without self consistency
        result = EM.SC_dilute(Km, Gm, Ki, Gi, f, 'stress')
        expected_result = (np.array([1.        , 1.18537201, 1.45510836, 1.88376754, 2.67045455]), np.array([1.        , 1.25214408, 1.67431193, 2.52595156, 5.14084507]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_SC_flex(self):
        ro=13/6
        Km, Gm=3*ro, 3 #
        Ki, Gi=Km*50, Gm*50  #  65, 30
        f= np.linspace(0,0.4,5)
        iter_n= 100
        result = EM.SC_flex(f,iter_n,Km,Ki,Gm,Gi)
        expected_result =(np.array([ 6.5       ,  7.72808652,  9.61307829, 12.88577277, 19.66902292]), np.array([ 3.        ,  3.7470363 ,  4.94410361,  7.10044198, 11.64810309]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_MT_average(self):
        f= np.linspace(0,1,5)
        Ki, Gi=5,10 #
        Km, Gm=37,45  #  
        result = EM.MT_average(f, Km, Gm,Km, Gm, Ki, Gi)
        expected_result = (np.array([ 5.        , 10.84269663, 17.83950617, 26.36986301, 37.        ]), np. array([10.        , 15.77381712, 23.02430955, 32.40103909, 45.        ]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
