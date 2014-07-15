#
# Author: Daniel G. A. Smith
# Created: 6/15/14
# Original content from:
# http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming
#

import numpy as np
from scipy import linalg as SLA

#Setup a few constats for the HF computation
nuclear =  8.002367061810450
ndocc = 5

print 'Reading in integrals...'

# Read in integrals
So = np.genfromtxt('S.dat', delimiter=',')
To = np.genfromtxt('T.dat', delimiter=',')
Vo = np.genfromtxt('V.dat', delimiter=',')
Io = np.genfromtxt('eri.dat', delimiter=',')

def make_array(arr):
    I1 = arr[:,0].astype(np.int) - 1
    I2 = arr[:,1].astype(np.int) - 1
    out = np.zeros((np.max(I1) + 1, np.max(I2) + 1))

    # 2 fold symmetry
    # Use numpy advanced indexing
    out[(I2,I1)] = arr[:,2]
    out[(I1,I2)] = arr[:,2]
    return np.matrix(out)

#### Normal integrals
S = make_array(So)
T = make_array(To)
V = make_array(Vo)

### ERI
sh = []
for x in range(4):
    sh.append(Io[:,x].astype(np.int) - 1)

### 8 fold symmetry
I = np.zeros(tuple(np.max(x)+1 for x in sh))
I[(sh[0],sh[1],sh[2],sh[3])] = Io[:,-1]
I[(sh[0],sh[1],sh[3],sh[2])] = Io[:,-1]

I[(sh[1],sh[0],sh[2],sh[3])] = Io[:,-1]
I[(sh[1],sh[0],sh[3],sh[2])] = Io[:,-1]

I[(sh[3],sh[2],sh[1],sh[0])] = Io[:,-1]
I[(sh[3],sh[2],sh[0],sh[1])] = Io[:,-1]

I[(sh[2],sh[3],sh[1],sh[0])] = Io[:,-1]
I[(sh[2],sh[3],sh[0],sh[1])] = Io[:,-1]

print '..Finished reading in integrals.\n'


# Compute Hcore 
H = T + V

# Orthogonalizer A = S^-1/2
A = np.matrix(SLA.sqrtm(S)).I
A = A.real

# Calculate initial core guess
# Using the matrix class
# * is equivalent to matrix multiplication
Hp = A * H * A
e,C2 = SLA.eigh(Hp)
C = A * C2 
D = C[:,:ndocc] * C[:,:ndocc].T

print('\nStarting SCF iterationations\n')
E    = 0.0
Enuc = nuclear
Eold = 0.0
maxiteration = 30
E_conv = 1E-10


for iteration in range(1, maxiteration + 1):

    # Fock Build
    J = np.einsum('pqrs,rs', I, D) 
    K = np.einsum('pqrs,qs', I, D)
    F = H + J * 2 - K

    E = np.einsum('ij,ij->', F + H, D) + Enuc

    # Roothaan Update
    print('@RHF Iteration %3d: Energy = %24.16f dE = %11.3E' % (iteration, E, E - Eold))
    if (abs(E - Eold) < E_conv):
        break

    Eold = E

    # New guess
    Fp = A * F * A
    e, C2 = SLA.eigh(Fp)
    C = A * C2
    D  = C[:,:ndocc] * C[:,:ndocc].T


print 'SCF Final Energy %5.10f' % E



