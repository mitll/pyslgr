# Get F0 library makefile

SRC=dp_f0.cpp f0.h f0_structs.h get_cands.cpp get_f0_fn.cpp sigproc.cpp
OBJ=$(SRC:.cpp=.o)

all: get_f0.a

get_f0.a : $(OBJ)
	ar rcv get_f0.a $(OBJ)

.cpp.o :
	g++ -O3 -fPIC -I../slgr_engine/ -c -std=c++0x -o $@ $<

clean :
	rm -f *.o *.a
