#!/usr/bin/env python3
# flake8: noqa
"""Setup script for the PyCIGAR repository."""
from os.path import dirname, realpath
from setuptools import find_packages, setup, Distribution
import setuptools.command.build_ext as _build_ext
import subprocess
from pycigar import __version__


def _read_requirements_file():
    """Return the elements in requirements.txt."""
    req_file_path = '%s/requirements.txt' % dirname(realpath(__file__))
    with open(req_file_path) as f:
        return [line.strip() for line in f]


"""class build_ext(_build_ext.build_ext):
    def run(self):
        try:
        except ImportError:
            subprocess.check_call([]))


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True
"""

setup(
    name='nr3',
    version=__version__,
    #distclass=BinaryDistribution,
    #cmdclass={"build_ext": build_ext},
    packages=find_packages(),
    description=("Powerflow with Newton-Raphson method"),
    long_description=open("README.md").read(),
    url="https://github.com/lbnl-cybersecurity/ceds-cigar",
    keywords=("distributed grid"
              "reinforcement-learning deep-learning python"),
    install_requires=_read_requirements_file(),
    zip_safe=False,
    scripts=['nr3/__init__.py',]
)
