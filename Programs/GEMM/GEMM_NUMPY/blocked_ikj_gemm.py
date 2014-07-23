#
# Author: Daniel Smith
# Date: 7/14/14
# Extremely simple numpy example of blocked dgemm code
# Will by slower than pure numpy, but this is much easier to see the blocking
#

import numpy as np

n = 30
bsize = 7

A = np.random.rand(n, n)
B = np.random.rand(n, n)

# All three dimensions of the block are equal size
num_blocks = n // bsize + 1

C = np.zeros((n, n))
for i in range(num_blocks):
    i *= bsize
    si = slice(i, i + bsize)
    for k in range(num_blocks):
        k *= bsize
        sk = slice(k, k + bsize)

        tmp_A = A[si, sk]  #Load tmp_A only once per loop
        for j in range(num_blocks):
            j *= bsize
            sj = slice(j, j + bsize)
              
            C[si, sj] += np.dot(tmp_A, B[sk, sj])


bench = np.dot(A,B)

works = np.allclose(C, bench)

if works:
    print 'The blocked DGEMM works!'
else:
    print 'The blocked DGEMM does not work.'
