###pySLGR 0.7.0

Preliminary release:
* This release has initial functionality and basic algorithms for running SLGR systems 
* Future releases will have hyper-parameter training, documentation, and example systems

Installation:
* Download the files from the git release tab 
* Type 'conda install pyslgr-0.7.0-py27_0.tar.bz2' to install pySLGR
* To run the examples untar the examples and models in the same directory
* 'cd examples' and then './example_signal.py' 

Building Cython interface only:
* Type 'make' to build
* Type 'make clean' to clean directories
* Type 'make eg' to see some examples -- more examples are available in the examples subdirectory

Building an Anaconda package:
* cd pyslgr
* conda build . 
* conda install pyslgr --use-local

