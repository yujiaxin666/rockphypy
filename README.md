[![Documentation Status](https://readthedocs.org/projects/rockphypy/badge/?version=latest)](http://rockphypy.readthedocs.io/en/latest/?badge=latest)[![PyPI version](https://badge.fury.io/py/rockphypy.svg)](https://badge.fury.io/py/rockphypy)[![Downloads](https://static.pepy.tech/badge/rockphypy)](https://pepy.tech/project/rockphypy)
---------


# [rockphypy](https://rockphypy.readthedocs.io/en/latest/)

If you have installed rockphypy before, upgrade to the latest version with:

```bash
pip install -U rockphypy
```

# Release Note rockphypy 0.0.3

Bug fixes
- `AVO.zoeppritz`: wrong results for nonzero angles (degree/radian bug and wrong medium in one term)
- `EM.PQ`, `EM.PQ_vectorize`: wrong P and Q for prolate inclusions (alpha > 1)
- `EM.PQ_vectorize`: no longer changes the input alpha array
- `EM.hudson`, `hudson_rand`, `hudson_ortho`, `hudson_cone`: wrong M term for cracks filled with a solid
- `EM.hudson_rand`: shear modulus correction was 15 times too large
- `EM.DEM`, `EM.Berryman_DEM`: corrected parameter order and inclusion shear modulus
- `Anisotropy.Bond_trans`: missing factor 2 in the Bond matrix
- `Empirical.han`: clay coefficient for Vp is 2.18 (Han et al., 1986, 40 MPa)
- `Fluid.White_Dutta_Ode`: numerically stable at very low and very high frequencies
- `QI.estimate_cem`: works with pandas 3

Compatibility
- Works with NumPy 2: requires KDE-diffusion >= 1.0.5
NumPy 2.0: older versions of rockphypy could fail to import with NumPy ≥ 2.0 because of the dependency KDE-diffusion 1.0.4. From rockphypy 0.0.3 we require KDE-diffusion ≥ 1.0.5, which works with NumPy 2. 

Documentation
- New documentation theme

# About

**rockphypy** is a Python-based open-source library dedicated to rock physics modeling. It caters to individuals interested in rock physics and empowers them to apply rock physics tools in their work and research effectively.

Built upon the foundations of Stanford SRB tools, rockphypy has successfully migrated and optimized a wide range of rock physics models from Stanford SRB Matlab tools to Python. Additionally, it extends the functionality by introducing new functions and practical workflows. Check the documentation and tutorials [here](https://rockphypy.readthedocs.io/).

## Installation

rockphypy is available through [PyPI] and may be installed using `pip`:

```text
pip install rockphypy
```

# Citation 

Want to cite rockphypy in your work? 

Bibtex: 

@article{YU2023101567, <br>
title = {rockphypy: An extensive Python library for rock physics modeling},<br>
journal = {SoftwareX},<br>
volume = {24},<br>
pages = {101567},<br>
year = {2023},<br>
issn = {2352-7110},<br>
doi = {https://doi.org/10.1016/j.softx.2023.101567},<br>
url = {https://www.sciencedirect.com/science/article/pii/S2352711023002637},<br>
author = {Jiaxin Yu and Tapan Mukerji and Per Avseth}}

### Timeline of development

![2timeline](https://github.com/yujiaxin666/rockphypy/assets/45630390/5689968c-7683-41e4-864a-0dca791a38a0)

