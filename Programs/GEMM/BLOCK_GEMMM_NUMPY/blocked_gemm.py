#
# Author: Daniel Smith
# Date: 7/14/14
# Extremely simple numpy example of blocked dgemm code
#

import numpy as np

n = 30
bsize = 10

A = np.random.rand(n, n)
B = np.random.rand(n, n)


num_blocks = n // bsize + 1

C = np.zeros((n, n))
for i in range(num_blocks):
    i *= bsize
    si = slice(i, i + bsize)
    for j in range(num_blocks):
        j *= bsize
        sj = slice(j, j + bsize)
        for k in range(num_blocks):
            k *= bsize
            sk = slice(k, k + bsize)
              
            # Different options (Pick one)
            C[si, sj] += np.dot(A[si, sk], B[sk, sj])
            # C[si, sj] += np.einsum('ij,jk', A[si, sk], B[sk, sj])
            # C[si, sj] += np.sum(A[si, sk][..., None] * B[sk, sj], axis=1)


bench = np.dot(A,B)

print np.allclose(C, bench)

