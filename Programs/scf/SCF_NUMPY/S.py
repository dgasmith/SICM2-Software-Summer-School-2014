import numpy as np
import pandas as pd
from scipy import linalg as SLA
import time
np.set_printoptions(precision=4, suppress=True)

nuclear =  8.002367061810450
ndocc = 5

So = pd.read_csv('S.dat',header=None).values 
To = pd.read_csv('T.dat',header=None).values
Vo = pd.read_csv('V.dat',header=None).values 
#Io = pd.read_csv('eri.dat',header=None).values

def make_array(arr):
    I1 = arr[:,0].astype(np.int)-1
    I2 = arr[:,1].astype(np.int)-1
    out = np.zeros((np.max(I1)+1,np.max(I2)+1))
    out[(I2,I1)] = arr[:,2]
    out[(I1,I2)] = arr[:,2]
#    out[np.tril_indices_from(out,k=-1)] = out[np.triu_indices_from(out,k=1)]
#    out = out + out.T - np.diag(out.diagonal())
    return out

#### Normal integrals
S = make_array(So)
T = make_array(To)
V = make_array(Vo)
#print S

#exit()
#Orthogonalizer A
A = np.matrix(SLA.sqrtm(S)).I
A = A.real

print A

e,C2 = SLA.eigh(S)
#C2 = C2

print e
e = e**(-0.5)
print e
e = np.diag(e)

print C2.dot(e).dot(C2.T)

