# This is in the public domain.
# Created by Daniel Smith 7/10/14
 
import glob
import pandas as pd
 
infiles = glob.glob('Data/*.out')
 
def read(infile):
    """Read 'infile' into a list where each element is a line of the file"""
    return open(infile, 'r').readlines()
 
def find(data, string, position=False, dtype=float):
    """Looks for the pattern 'string' in data"""
 
    fdata = filter(lambda x: string in x,data)
 
    if position is False:
        return list(fdata)
    else:
        return [dtype(x.split()[position]) for x in fdata]
 
 
output = []
for inp in infiles:
 
    #Grab systematic indices from the filename
    name = inp.split('/')[-1].split('+')
    index = name[1].split('_')
 
    #Grab RKS energies
    data = read(inp)
    energies = find(data, '!RKS STATE 1.1 Energy', -1)
    energy = (energies[0] - energies[1] - energies[2])*627.509
 
    #append to output
    output.append(index + [energy])
 
#Create and print pandas DataFrame
output = pd.DataFrame(output)
output.columns = ['r','d','theta','x','y','z','BLYP/QZVP']
 
print output
