from opentps.core.processing.doseCalculation.photons._utils import shiftBeamlets2,shiftBeamlets 
import os
import numpy as np
import pickle
import scipy.sparse as sp
import copy

def loadBeamlets(file_path):
    from opentps.core.data._sparseBeamlets import SparseBeamlets
    return loadData(file_path, SparseBeamlets)

def loadData(file_path, cls):
    with open(file_path, 'rb') as fid:
        tmp = pickle.load(fid)
    data = cls()
    data.__dict__.update(tmp)
    return data

angles = [0]*int(3821/3) + [np.pi/2]*int(3821/3)
angles+= [np.pi/4]*(3821-len(angles))
doseGridSize = [100,100,100]

sp1 = sp.csc_matrix(([1.1, 2.2, 3.3], ([1000,1001,1002], np.zeros(3))), shape=(np.prod(doseGridSize), 1),dtype=np.float32)   
sp2 = sp.csc_matrix(([5.1, 6.2, 7.3], ([1000,1001,1002], np.zeros(3))), shape=(np.prod(doseGridSize), 1),dtype=np.float32)   
sparseBeamlets = sp.hstack([sp1, sp2])

# shift = [4.1, -5.1, -4.2]
shift = [1.4, -1.2 , -2.3]
BeamletMatrix = shiftBeamlets2(sparseBeamlets, doseGridSize, shift, [2.3, 1.3]) ### 
shift = [1.4, -1.2 , -2.3]
# BeamletMatrix1 = shiftBeamlets(sparseBeamlets, doseGridSize, shift, [2.3, 1.3]) ### 
print(BeamletMatrix.nonzero())
print(BeamletMatrix[BeamletMatrix.nonzero()])

a=0