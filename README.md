# pySLGR 1.0

A Python tool for Speaker, Language, and Gender Recognition (SLGR) and general meta-data extraction.

## Preliminary release:
* This release has initial functionality and basic algorithms for running SLGR systems 
* Future releases will have hyper-parameter training, documentation, and example systems

## Requirements:
* Linux (Ubuntu)
* Boost libraries
* Anaconda v4.3.1 for Python 2.7.x

## Documentation:
[Online documentation](https://mitll.github.io/pyslgr) is avaialable at: [https://mitll.github.io/pyslgr](https://mitll.github.io/pyslgr)

## Installation:
* Download the files from the git release tab 
* Type 'conda install pyslgr-0.7.2-py27_0.tar.bz2' to install pySLGR
* To run the examples untar the examples and models in the same directory
* 'cd examples' and then run examples with the command 
		python <example_name.py>
	For instance to run the example_signal example type:
		python example_signal.py

## Building pySLGR from Source Code
### 1. Build Cython Interface:
* Type 'make' to build
* Type 'make clean' to clean directories
* Type 'make eg' to see some examples -- more examples are available in the examples subdirectory

### 2. Build Anaconda package:
* Install conda-build package (version 2.1.10) in Anaconda if not installed already. Use the following command:
		conda install conda-build=2.1.10 
* cd <installation_dir>/pyslgr/pyslgr
* Run the following command:
      	      conda build . 

## Installing Conda Package After Build
* Use the following command:
		conda install pyslgr --use-local
## Citation
Please use this DOI number reference, published on Zenodo, when citing the software:

[![DOI](https://zenodo.org/badge/76076310.svg)](https://zenodo.org/badge/latestdoi/76076310)

## Disclaimer

DISTRIBUTION STATEMENT A. Approved for public release: distribution unlimited.

© 2020 MASSACHUSETTS INSTITUTE OF TECHNOLOGY

    Subject to FAR 52.227-11 – Patent Rights – Ownership by the Contractor (May 2014)
    SPDX-License-Identifier: MIT

This material is based upon work supported by the US Air Force under Air Force Contract No. FA8721-05-C-0002 and/or FA8702-15-D-0001. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the US Air Force.

The software/firmware is provided to you on an As-Is basis

