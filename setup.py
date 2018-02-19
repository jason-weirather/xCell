from setuptools import setup, find_packages
from codecs import open
from os import path
import sys

if sys.version_info < (3,4) and sys.version_info < (2.7):
   sys.exit("Error: You are using Python "+str(sys.version_info)+"; Please use Python  3.4 or 2.7 or better\n")

this_folder = path.abspath(path.dirname(__file__))
with open(path.join(this_folder,'README.md'),encoding='utf-8') as inf:
  long_description = inf.read()

setup(
  name='xcell',
  version='1.1.0.0',
  description='Python CLI and module for running the xCell package with Python Pandas inputs and outputs.',
  long_description=long_description,
  url='https://github.com/jason-weirather/xCell',
  author='Jason L Weirather',
  author_email='jason.weirather@gmail.com',
  license='Apache License, Version 2.0',
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: Apache Software License'
  ],
  keywords='bioinformatics, R, enrichment, xCell, ssGSEA, microenvironment, RNAseq, immune infiltrates',
  packages=['xcell'
           ],
  package_data={'xcell':['*.r']},
  install_requires=['pandas'],
  entry_points = {
    'console_scripts':['xcell=xcell:__cli']
  }
)

