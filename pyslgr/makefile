# pyslgr make for Linux
# 

PLIB_EXT=so
LIBS=../slgr_engine/slgr_lib.a ../get_f0_lib/get_f0.a 

all : LLSignal.$(PLIB_EXT) LLFeatures.$(PLIB_EXT) MFCCFeatures.$(PLIB_EXT) GMMModel.$(PLIB_EXT) 

LLSignal.$(PLIB_EXT) : LLSignal.pyx Signal.pxd LLSignal.pxd $(LIBS)
	rm -f pyslgr/LLSignal.$(PLIB_EXT)
	python setup_signal.py build_ext -i

LLFeatures.$(PLIB_EXT) : LLFeatures.pyx Features.pxd $(LIBS)
	rm -f LLFeatures.$(PLIB_EXT)
	python setup_features.py build_ext -i

MFCCFeatures.$(PLIB_EXT) : MFCCFeatures.pyx LLFeatures.pyx LLFeatures.pxd MFCCFeatures.pxd MFCC_Features.pxd $(LIBS)
	rm -f MFCCFeatures.$(PLIB_EXT)
	python setup_mfcc.py build_ext -i

GMMModel.$(PLIB_EXT) : GMMModel.pyx GMMModel.pxd MFCCFeatures.$(PLIB_EXT) LLFeatures.$(PLIB_EXT) $(LIBS)
	rm -f GMMModel.$(PLIB_EXT)
	python setup_gmm.py build_ext -i

clean: 
	\rm -rf build
	\rm -f *.cpp
	\rm -f *.so *.pyc

