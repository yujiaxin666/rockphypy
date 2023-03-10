# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding

from os import path
import re
# The directory containing this file
HERE = path.abspath(path.dirname(__file__))



def read(fname):
    return open(path.join(HERE, fname)).read()


# Get version
with open(path.join(HERE, "rockphypy/__init__.py"), encoding="utf-8") as file:
    VERSION = re.search(r"__version__ = \"(.*?)\"", file.read()).group(1)

# Get the long description from the README file
# with open(path.join(HERE, 'README.md'), encoding='utf-8') as f:
#     long_description = f.read()

# This call to setup() does all the work
setup(
    name="rockphypy",
    version=VERSION,
    description="Python rock physics toolbox  ",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    #url="https://pypi.org/project/RockPhyRoll/0.0.1/",
    author="Jiaxin Yu",
    author_email="yujiaxin666@outlook.com",
    license="GNU GPLv3",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",

    ],
    packages=["rockphypy"],
    include_package_data=True,
    install_requires=["numpy>=1.14.2", "scipy>=1.0.1", "matplotlib>=2.2.0"]
)