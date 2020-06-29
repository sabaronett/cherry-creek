import numpy as np

# load memory output data
data = np.loadtxt('output/seqtimes.txt') # return (N, 2) array
loop = data[:, 0]                        # return only 1st col
time = data[:, 1]                        # return only 2nd col

print('---Time Performance Statistics---')
print('Total runs:   %d'%(loop.size))
print('   Average:   %f s'%(np.average(time)))
print('     Error: ± %f s'%(np.std(time)))
