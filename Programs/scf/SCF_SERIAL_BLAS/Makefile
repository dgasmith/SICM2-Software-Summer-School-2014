#CXX = g++
#CXXLIBS = -L/home/guest40/lib -llapack3 -lblas3 -lgfortran
CXX = icpc
CXXLIBS = -mkl
CXXFLAGS = -g -Wall -std=c++11

scf: scf.o diag.o mmult.o
	$(CXX) $^ $(CXXLIBS) -o $@

clean:
	rm -f *.o scf
