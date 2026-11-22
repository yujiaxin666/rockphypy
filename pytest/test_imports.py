# Check that every module of rockphypy can be imported.
# Problems at import time (for example a new numpy, pandas or matplotlib version
# rejecting something that runs when a module is loaded) make the whole package unusable.
import importlib
import pkgutil

import pytest
import rockphypy

MODULES = sorted(m.name for m in pkgutil.iter_modules(rockphypy.__path__, "rockphypy."))


@pytest.mark.parametrize("name", MODULES)
def test_import_module(name):
    importlib.import_module(name)


@pytest.mark.parametrize("name", ["AVO", "Anisotropy", "BW", "EM", "Empirical",
                                  "Fluid", "GM", "Permeability", "QI", "utils"])
def test_public_names(name):
    # the names users import with: from rockphypy import ...
    assert hasattr(rockphypy, name)