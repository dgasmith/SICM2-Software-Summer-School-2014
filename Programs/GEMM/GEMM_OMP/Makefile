CXX      = icpc
CXXFLAGS = -mkl -std=c++11 -O3 -xHOST -openmp -restrict
CPPFLAGS = -DHAVE_MKL
LDFLAGS  = -mkl

# uncomment if want to use Eigen
# EIGENPATH = /gpfs/software/eigen-eigen-6b38706d90a9
EIGENPATH = /Users/daniel/Downloads/eigen-eigen-6b38706d90a9
CPPFLAGS  += -I$(EIGENPATH) -DHAVE_EIGEN

default:: gemm

# gemm:: gemm.o Kernals/standard_gemm.o
gemm:: gemm.o Kernals/standard_gemm.o Kernals/blocked_gemm.o
	$(CXX) $(LDFLAGS) $^ -o $@

clean:
	rm -f *.o Kernals/*.o gemm

