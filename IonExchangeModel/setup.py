# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 06:46:48 2020

@author: LHaupert
"""

"""setup.py: setuptools control."""


from setuptools import setup

with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")


setup(
    name = "ion-exchange-models",
    packages = ["ixpy",],
    entry_points = {
        "console_scripts": ['ixpy = ixpy.__main__:main']
        },
    version = "0.1.0",
    description = "Ion Exchange Kinetic Models",
    long_description = long_descr,
    author = "Levi Haupert, Boris Datsov, Jonathan Burkhardt",
    author_email = "haupert.levi@epa.gov",
    license='MIT',
    url='https://github.com/USEPA/Water_Treatment_Models',
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering"]
    )