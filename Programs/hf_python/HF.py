#############################################################
#                 Wallace Derricotte                        #                
#                 Hartree-Fock Theory                       #
#                Python Coding Example                      #
#############################################################


import subprocess
import os
from numpy import *
from numpy.linalg import *

# Define function that reads one-electron integrals
def read_one_electron_integrals(filename):
    ints_lines = open(filename,"r").readlines()
    n = len(ints_lines[0].split())
    ints = zeros( (n,n) )
    for i,line in enumerate(ints_lines):
        for j,value in enumerate(line.split()):
            ints[i][j] = float(value)
    return ints

# Define function that reads two-electron integrals
def read_two_electron_integrals(filename,nao):
    ints_lines = open(filename,"r").readlines()
    ints = zeros( (nao,nao,nao,nao) )
    for line in ints_lines:
        line_split = line.split()
        p,q,r,s = [int(t) for t in line_split[0:4]]
        value = float(line_split[-1])
        ints[p][q][r][s] = ints[q][p][r][s] = value
        ints[p][q][s][r] = ints[q][p][s][r] = value
        ints[r][s][p][q] = ints[r][s][q][p] = value
        ints[s][r][p][q] = ints[s][r][q][p] = value
    return ints

##############################################################
#                  User Input Block                          #
##############################################################
nuc_nuc = "vnn" #File containing Nuc-Nuc Repulsion
One_e = "one-electron" #File Containing 1e- integrals
Two_e = "two-electron" #File Containing 2e- integrals
S = "overlap" #File Containing overlap integrals
e_tot = 10 #Specify total number of electrons in your molecule
##############################################################


# Open files to read Nuc-Nuc repulsion value, also opens output
# file for writing.
o = open(nuc_nuc,"r")
f = open("output.dat","w")
vnnstr = o.readline()
vnn = float(vnnstr)
f.write("Nuclear-Nuclear Repulsion: %f \n \n" %vnn)


# Read one-electron integrals from file and store in Matrix
OneElec = read_one_electron_integrals(One_e)
f.write("Core Hamiltonian \n %s \n" %OneElec)

# Number of Atomic Orbitals can be obtained from the
# dimensions of the Core Hamiltonian(OneElec)
nao = len(OneElec)

# Read two-electron integrals from file and store in Matrix
TwoElec = read_two_electron_integrals(Two_e, nao)

# Check integrals by computing sum of (pq|rs)^2 and abs(pq|rs)
a = 0
for i in range(nao):
    for j in range(nao):
        for k in range(nao):
            for l in range(nao):
                a += (TwoElec[i][j][k][l])**2
f.write("\n \nSum of (pq|rs)^2: %f \n" %a)
b = 0
for i in range(nao):
    for j in range(nao):
        for k in range(nao):
            for l in range(nao):
                b += abs(TwoElec[i][j][k][l])
f.write("Sum of abs(pq|rs): %f \n" %b)

# Read in overlap integrals S and store in Matrix
Overlap = ()
Overlap = genfromtxt(S, usecols = range(nao))
f.write("\n \nOverlap Matrix \n %s" %Overlap)

# Get Eigenvalues of S matrix
Eigen_Overlap = eigh(Overlap)
eigvecs = Eigen_Overlap[1]

# Form s^(-1/2) eigenvalue matrix
s_half = zeros( (nao,nao) )
S_half = zeros( (nao,nao) )
Left_Hand = zeros( (nao,nao) )
for i in range(nao):
    for j in range(nao):
        if (i==j):
            s_half[i][j] = 1.0/((Eigen_Overlap[0][i])**0.5)
        else:
            s_half[i][j] = 0.0
Left_Hand = dot(eigvecs,s_half)
S_half = dot(Left_Hand,eigvecs.T)
f.write("\n\nX Matrix \n %s \n" %S_half)


# A quick check to see if S_half is correct, Y should be
# the identity matrix.
Y = zeros( (nao,nao) )
L = zeros( (nao,nao) )
L = dot(S_half.T, Overlap)
Y = dot(L,S_half)
f.write("\nXSX Matrix \n %s" %Y)


# Set interation counter, Density Matrix, and initial energy
# all to zero. DensityOld will only be used as a placeholder
# for the density matrix from the previous iteration.
count = 0
Density = zeros( (nao,nao) )
DensityOld = zeros( (nao,nao) )
initial_energy = 0.0

# The energy change and RMS are set to absurdly high values
# so the loop can begin, and to easily spot early errors in
# the loop.
del_energy = 10000000
RMS = 10000000

# Begining the iterative procedure
while (abs(del_energy) > 1e-9 or abs(RMS) > 1e-9):
    f.write("\n\n*******************")
    f.write("\n Iteration: %d \n" %count)
    f.write("*******************\n")
    if count > 0:
        f.write("\n\n         HERE WE GO AGAIN!! Another spin 'round the wheel of convergence\n")
    else:
        f.write("\n\n                      Let's Begin!! \n")
# Increase count by one
    count = count+1
#Form G Matrix
    G = zeros( (nao,nao) )
    for i in range(nao):
        for j in range(nao):
            a = 0
            for k in range(nao):
                for l in range(nao):
                    a += Density[k][l]*(2.0*TwoElec[i][j][k][l]
                                        - TwoElec[i][k][j][l])
                    G[i][j] = a
    f.write("\nG Matrix \n %s" %G)
# Form Fock Matrix
    Fock = zeros( (nao,nao) )
    Fock = OneElec + G
    f.write("\n\nFock Matrix(F)\n %s \n" %Fock)
# Compute the right hand side of the energy expression
# dhf = D(H + F)
    dhf = 0
    for i in range(nao):
        for j in range(nao):
            dhf += Density[i][j]*(OneElec[i][j] + Fock[i][j])
# Compute Energy
    energy = vnn + dhf
    del_energy = initial_energy - energy
# Transform Fock Matrix to orthonormal AO basis
    FockPrime = zeros( (nao,nao) )
    L = zeros( (nao,nao) )
    L = dot(S_half,Fock)
    FockPrime = dot(L,S_half)
    f.write("\n\nF' Matrix \n %s \n" %FockPrime)
# Get eigenvalues/vectors of new Fock Matrix (CoeffPrime)
    Eigen_FockPrime = eigh(FockPrime)
    f.write("\n\nEigenvalues of F' \n %s" %Eigen_FockPrime[0])
    CoeffPrime = zeros( (nao,nao) )
    Coeff = zeros( (nao,nao) )
    CoeffPrime = Eigen_FockPrime[1]
    f.write("\n\nC Matrix \n %s" %CoeffPrime)
# Transform coefficients back to original AO basis
    Coeff = dot(S_half,CoeffPrime)
    f.write("\n\nC' Matrix \n %s" %Coeff)
    nmo = e_tot/2
# Form New Density Matrix
    for i in range(nao):
        for j in range(nao):
            a = 0
            for z in range(nmo):
                a += Coeff[i][z]*Coeff[j][z]
                Density[i][j] = a
    f.write("\n\nD Matrix' \n %s" %Density)
    rms = 0
# Compute RMS
    for i in range(nao):
        for j in range(nao):
            rms += (Density[i][j] - DensityOld[i][j])**2
    RMS = sqrt(rms)
# Display Results of the iteration
    f.write("\n \nSCF Energy:%.10f \n" %energy)
    f.write("Energy Change: %.10f\n" %del_energy)
    f.write("Root Mean Square Deviation: %.10f" %RMS)
    RMS = 0
    rms = 0
    initial_energy = energy
    for i in range(nao):
        for j in range(nao):
            DensityOld[i][j] = Density[i][j]
# Ensure code doesn't get stuck in senseless loop
    if count > 100:
        break

# Display Final Results
f.write("\n\n\n\n              **************************************************")
f.write("\n              ********  FINAL SCF ENERGY: %.10f *******" %initial_energy)
f.write("\n              **************************************************")
f.write("\n\n\n   'If you don't finish what you start, your success rate will always be zero.' ")

#DONE!!





