# simple makefile to compile pyslgr

all: 
	cd get_f0_lib ; make
	cd slgr_engine ; make
	cd pyslgr ; make

clean:
	cd get_f0_lib ; make clean
	cd slgr_engine ; make clean
	cd pyslgr ; make clean
	cd examples ; rm -rf tmp/

eg:
	cd examples ; \
	./example_signal.py ; \
	./example_mfcc.py ; \
	./example_mfcc_fb.py ; \
	./example_lid_gmm.py
