#
# Author: Daniel Smith
# Date: 7/14/14
# Extremely simple numpy example of j,k unrolled dgemm code
# Will by slower than pure numpy, but this is much easier to see
#

import numpy as np

n = 40

A = np.random.rand(n, n)
B = np.random.rand(n, n)
C = np.zeros((n, n))

for i in range(n):
    for j in range(0, n, 2):
        cij_1 = np.zeros(2)
        cij_2 = np.zeros(2)

        for k in range(0, n, 2):
            sk = slice(k, k+2)

            tmp_a = A[i, sk]     #Notice how we only have to load array A once!
            tmp_b1 = B[sk, j]
            tmp_b2 = B[sk, j+1]

            cij_1 += tmp_a * tmp_b1
            cij_2 += tmp_a * tmp_b2

        C[i, j] += cij_1[0] + cij_1[1]
        C[i, j+1] += cij_2[0] + cij_2[1]

bench = np.dot(A,B)

works = np.allclose(C, bench)

if works:
    print 'The blocked DGEMM works!'
else:
    print 'The blocked DGEMM does not work.'
