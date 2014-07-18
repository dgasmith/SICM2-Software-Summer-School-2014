## GEMM BLAS, Plain, and Blocked routines
`C[i][j] = A[i][k] * B[k][j]`

#### Blocked GEMM
The primary effort of this code went into the development of the blocked gemm routine. 
With the current setting `block8` provides the best performance with a ikj loop ordering
and with the kj loops manually unrolled 4 times with blocks of shape `4x8x1`.
Other algorithms are present in the code and can be called by simply managing the comments
in the inner loop of the `blocked_gemm` routine in `Kernals/blocked_gemm.cc`. 


#### loop_gemm.py
Provides an easy to use python interface for timing the gemm program. 

