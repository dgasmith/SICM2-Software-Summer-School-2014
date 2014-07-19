EIGENPATH = /Users/adrianmo/Documents/OSULocal/Research/code/SICM2-Software-Summer-School-2014/
LIBINT2PATH = /Users/adrianmo/Documents/OSULocal/Research/code/SICM2-Software-Summer-School-2014/libint-2.0.3
LIBINT2INCLUDES = -I$(LIBINT2PATH)/include -I$(LIBINT2PATH)/include/libint2

CXX = g++
CXXFLAGS = -O2 -march=native -std=c++11 -fopenmp
CPPFLAGS = -I$(EIGENPATH) $(LIBINT2INCLUDES)
LDFLAGS = -L$(LIBINT2PATH)/lib -lint2

default:: scf

scf: scf.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) scf.cc -o scf $(LDFLAGS)

clean:
	rm -f *.o scf
