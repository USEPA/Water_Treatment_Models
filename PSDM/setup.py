# -*- coding: utf-8 -*-
"""
Created on Mon May 24, 2021
@author: UCChEJBB
"""

"""setup.py: setuptools control."""


from setuptools import setup
import os

cwd = os.getcwd()
# print(f"ixpy @ {cwd}../../Water_Treatment_Models/IonExchangeModel/ixpy")

with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")


setup(
    name = "PSDM-models",
    packages = ["PSDM",],
    # entry_points = {
    #     "console_scripts": ['ixpy = ixpy.__main__:main']
    #     },
    version = "0.1.0",
    description = "Pore and Surface Diffusion Models",
    long_description = long_descr,
    author = "Jonathan Burkhardt, Levi Haupert, Boris Datsov",
    author_email = "burkhardt.jonathan@epa.gov",
    license='MIT',
    url='https://github.com/USEPA/Water_Treatment_Models',
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering"],
    install_requires = [
        "matplotlib",
        "numpy",
        "pandas",
        "scipy",
        "openpyxl",
        "xlsxwriter",
        f"ixpy @ ..\\..\\Water_Treatment_Models\\IonExchangeModel\\"   
        ## {os.getcwd()}\\
    ],
    )