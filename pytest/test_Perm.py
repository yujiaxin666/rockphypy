#%%
#from import_path import rockphypy
from rockphypy import Permeability
import pytest
import numpy as np

class TestPerm:
    def test_Kozeny_Carman(self):
        result = Permeability.Kozeny_Carman(0.3,0.0001)
        expected_result = 3.0612244897959185e-12
        assert np.allclose(result, expected_result, rtol=1e-6, atol=1e-6)