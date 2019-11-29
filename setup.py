#!/usr/bin/env python

from setuptools import setup

# Version number
major = 2019
minor = 1

setup(
    name="surfaise",
    version="{}.{}".format(major, minor),
    description="Surfaise: Solving PDEs on parametrized surfaces in FEniCS.",
    author="Gaute Linga and Bjarke Frost Nielsen",
    author_email="gaute.linga@mn.uio.no",
    url="https://github.com/gautelinga/surfaise.git",
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python ",
        "License :: OSI Approved :: MIT",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Software Development :: Libraries :: Python Modules"],
    packages=["surfaise",
              "surfaise.common",
              "surfaise.utilities",
              "surfaise.analysis_scripts"],
    package_dir={"surfaise": "surfaise"},
    entry_points={"console_scripts": [
        "surfaise-visualize=surfaise.run_surfaise:visualize",
        "surfaise-postprocess=surfaise.run_surfaise:postprocess"]},
)
