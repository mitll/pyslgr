#
# makefile for SLGR engine
#

all : slgr_lib.a 

SRC=cluster_seg.cpp feat_general.cpp fft.cpp gmm.cpp segmentation.cpp signal_general.cpp feat_mfcc.cpp \
	frame.cpp ipdf_expansion.cpp speech_tools.h gmmsad.cpp get_f0.cpp
OBJ=$(SRC:.cpp=.o)

slgr_lib.a : $(OBJ)
	ar rcv slgr_lib.a $(OBJ)

.cpp.o : 
	g++ -fPIC -O3 -c -std=c++0x -o $@ $<

clean :
	rm -f *.o *.a
