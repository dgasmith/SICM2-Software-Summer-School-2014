#
# Author: Daniel Smith
# Date: 7/14/14
# Extremely simple numpy example of blocked dgemm code
# Will by slower than pure numpy, but this is much easier to see the blocking
#

import numpy as np

n = 30
bsize = 10

A = np.random.rand(n, n)
B = np.random.rand(n, n)

# All three dimensions of the block are equal size
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
              
            # Different options (Pick one):

            C[si, sj] += np.dot(A[si, sk], B[sk, sj])
            # C[si, sj] += np.einsum('ij,jk', A[si, sk], B[sk, sj])
            # C[si, sj] += np.sum(A[si, sk][..., None] * B[sk, sj], axis=1)


bench = np.dot(A,B)

works = np.allclose(C, bench)

if works:
    print 'The blocked DGEMM works!'
else:
    print 'The blocked DGEMM does not work.'
