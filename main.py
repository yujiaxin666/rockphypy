#%%
from rockphypy import Anisotropy
from rockphypy import GM

# %%
from rockphypy import AVO
from rockphypy import BW
from rockphypy import utils
from rockphypy import Permeability


# %%
Permeability.Bernabe(0.1,0.001,10,10)
# %%
utils.crack_por(0.001,1)
# %%

GM.johnson(37,30,6,0.1,0.1,0.2,'together')


import numpy as np
# %%
Anisotropy.Backus( [0.1,0.9],[2,8],[20,30] )
# %%
