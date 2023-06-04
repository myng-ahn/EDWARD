import os
from setuptools import setup, find_packages

# version-keeping code based on pybedtools
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 1
MIN = 0
REV = 1
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'edward/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )


setup(
    name='edward',
    version=VERSION,
    description='CSE 185 Demo Project',
    author='MJM',
    author_email='ADD',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "edward=edward.edward:main"
        ]
    }
)
