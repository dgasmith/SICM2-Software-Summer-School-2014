CXX      = icpc
CXXFLAGS = -std=c++11 -O3 -xHOST
LDFLAGS  = -mkl

default:: integrate

integrate:: integrate.o integrate_kernals.o
	$(CXX) $(LDFLAGS) $^ -o $@

clean:
	rm -f *.o integrate
