#comment out these three lines to build regular version. Uncomment to build MPI
PROGRAM = mmmpi
PROGSRC = mmmpi.cc
CXX = mpiicpc

#comment out these three lines to build MPI version. Uncomment to build Regular.
PROGRAM = mm
PROGSRC = mm.cc
CXX = icpc

CXXFLAGS = -g -std=c++0x -O3 -xHOST -openmp 
CXXLIBS = -mkl 
OBJECTS = $(PROGSRC:%.cc=%.o)

$(PROGRAM): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ $(CXXLIBS) -o $@

%.o: $.cc
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(PROGRAM)
