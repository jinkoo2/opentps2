from opentps.core.processing.doseCalculation.photons._utils import shiftBeamlets2,shiftBeamlets 
import os
import numpy as np
import pickle
import time


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
file = loadBeamlets('/home/luciano/ResearchData/DataSets/ICRU-Prostate/output/ProbabilisticProstate/SparseMatrices_Simulation/SM_nominal.pkl')
shift = [1.4, -1.2 , -2.3]
t0 = time.time()
BeamletMatrix_1 = shiftBeamlets2(file._sparseBeamlets, file.doseGridSize, shift, angles) ### 
# shift = [1.4, -1.2 , -2.3]
# BeamletMatrix = shiftBeamlets(file._sparseBeamlets, file.doseGridSize, shift, angles) ### 
print(time.time()-t0)
# a=0
# print(BeamletMatrix_1)