#from import_path import rockphypy
from rockphypy import Empirical

import numpy as np

class TestEmp:
    def test_krief(self):
        K0, G0= 47,44
        phi= np.linspace(1e-6,0.4,5 )
        result= Empirical.krief(phi, K0,G0)
        expected_result= (np.array([17.66661367, 12.4344996 ,  7.65139849,  3.83078846,  1.37376   ]), np.array([43.999868  , 30.96894241, 19.05631323,  9.54083163,  3.42144   ]), np.array([46.999859  , 33.08046121, 20.35560731, 10.19134288,  3.65472   ]))
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_esti_VS(self):
        Vp= np.array([3000,3500,2300])
        vsh= np.array([0.1,0.2,0.15])
        result = Empirical.esti_VS(Vp, vsh)
        expected_result= np.array([1544.70332922, 1931.5034334 ,  979.50251832])
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)
    def test_Cp(self):
        assert np.allclose(Empirical.Cp(0.4), 8.639999999999999, rtol=1e-6, atol=1e-6)

    