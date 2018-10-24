#!/usr/bin/env python
#coding=utf-8

import os
from setuptools import setup

setup(
    name = "ephpatic coupling library",
    version = "0.0.1",
    author = "Bartosz Telenczuk",
    author_email = "bartosz.telenczuk@gmail.com",
    description = ("model ephaptic interactions"),
    license = "HBP",
    keywords = "electric field, LFP, visual cortex, ephaptic, computational neuroscience",
    url = "http://github.com/btel",
    packages=['ei_model'],
    install_requires=[
        'attrs>=17.3.0',
        'ipykernel>=4.7.0',
        'ipython>=6.1.0',
        'jupyter-client>=5.1.0',
        'jupyter-console>=5.2.0',
        'jupyter-core>=4.3.0',
        'matplotlib==2.1.0',
        'nbconvert>=5.3.1',
        'nbformat>=4.4.0',
        'notebook>=5.2.2',
        'numpy==1.13.3',
        'scipy==0.19.1',
        'seaborn==0.8',
        'svgutils==0.3.1'   
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities"
    ],
)
