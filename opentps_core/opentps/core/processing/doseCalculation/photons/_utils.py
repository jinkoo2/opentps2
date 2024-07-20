import numpy as np
import opentps.core.io.CCCdoseEngineIO as CCCdoseEngineIO
import scipy.sparse as sp
from scipy import ndimage
import pycuda.autoinit
import pycuda.gpuarray as gpuarray
from pycuda.compiler import SourceModule
from pycuda.autoinit import context
import pycuda.driver as cuda
from scipy.ndimage import shift, gaussian_filter
import scipy.sparse as sp
from scipy.sparse import csc_matrix
from opentps.core.data.images import DoseImage
from opentps.core.data.plan._photonPlan import PhotonPlan
import ctypes
import os
import psutil

def getConvolveNonZeroElements(kernel_size, nonZeroIndexes, image_size):
    nonZeroIndexes_convolved = []
    for i in range(kernel_size):
        for j in range(kernel_size):
            for k in range(kernel_size):
                nonZeroIndexes_convolved.append(nonZeroIndexes + CCCdoseEngineIO.convertTo1DcoordFortran([i - kernel_size//2, j - kernel_size//2, k - kernel_size//2], image_size))
    return np.unique(nonZeroIndexes_convolved)

def convolveVoxel(sparse,index,kernel, image_size):    
    kernel_size = kernel.shape[0]
    convolution = 0
    for i in range(kernel_size):
        for j in range(kernel_size):
            for k in range(kernel_size):
                index_shifted = index + CCCdoseEngineIO.convertTo1DcoordFortran([i - kernel_size//2, j - kernel_size//2, k - kernel_size//2], image_size)
                convolution += sparse[index_shifted,0] * kernel[i,j,k]
    return convolution

mod=SourceModule("""
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

__device__ float getElement(int *index, int size, int target, float *sparse);
__device__ float linearInterpolation(float x, float xi0, float xi1, float xi2, float yi0, float yi1, float yi2);           
__device__ int sign(float x);

__global__ void shiftSparse(float *sparse, float *sparseShifted, int *index, int *indexShifted, int *shift, int *gridSize, float shiftValue, int indexNumberOfElements) {
    int start = threadIdx.x + blockIdx.x*blockDim.x;
    if (start < indexNumberOfElements) {
        indexShifted[2 * start] =  index[start] + shift[0];
        if (shiftValue != 0.0) {
            int index0 = index[start] + shift[1];
            float value0 = getElement(index, indexNumberOfElements, index0, sparse);

            float value1 = sparse[start];
                 
            int index2 = index[start] + shift[2];
            float value2 = getElement(index, indexNumberOfElements, index2, sparse);
            
            sparseShifted[2 * start] = linearInterpolation(shiftValue * -1, -1, 0, 1, value0, value1, value2);
            sparseShifted[2 * start + 1] = linearInterpolation(sign(shiftValue) - shiftValue, -1, 0, 1, value0, value1, value2);
                 
            if (value2 != 0.0)
              sparseShifted[2 * start + 1] /= 2;
            if (value0 != 0.0)
              sparseShifted[2 * start] /= 2;
            
            indexShifted[2 * start + 1] = indexShifted[2 * start] + shift[2];
      }
      else {
          sparseShifted[2 * start] = sparse[start];
      }          
    }         
}
__device__ float getElement(int index[], int size, int target, float sparse[]) {
    int left = 0;
    int right = size - 1;
    while (left <= right) {
        int mid = left + (int) ((right - left) / 2);
        if (index[mid] == target) {
            return sparse[mid]; // Element found, return its index
        } else if (index[mid] < target) {
            left = mid + 1; // Search the right half
        } else {
            right = mid - 1; // Search the left half
        }
    }
    return 0; // Element not found
}        

__device__ float linearInterpolation(float x, float xi0, float xi1, float xi2, float yi0, float yi1, float yi2) {
              
    if ((xi0 <= x) && (x <= xi1)) {
        return (yi1 - yi0) / (xi1 - xi0) * (x - xi0) + yi0; 
    } else if ((xi1 < x) && (x <= xi2)) {
        return (yi2 - yi1) / (xi2 - xi1) * (x - xi1) + yi1; 
    }      
    else {
        return -1000;
    }
}          

__device__ int sign(float x) {
    if (x>0) {
        return 1; 
    } else if (x<0) {
        return -1; 
    }      
    else {
        return 0;
    }
}   
""")

def correctShift(setup, angle):
    return np.array([(setup[0] * np.cos(angle) - setup[1] * np.sin(angle)) * np.cos(angle), (setup[0] * np.cos(angle) - setup[1] * np.sin(angle)) * np.sin(angle), setup[2]])

def shiftBeamlets(sparseBeamlets, gridSize,  scenarioShift_voxel, beamletAngles_rad):
    scenarioShift_voxel[2]*=-1 ### To have the setup error in LPS. Check because some signs problem
    scenarioShift_voxel[1]*=-1 ### To have the setup error in LPS. Check because some signs problem
    nbOfBeamlets = sparseBeamlets.shape[1]
    nbOfVoxelInImage = sparseBeamlets.shape[0]
    assert(len(beamletAngles_rad), nbOfBeamlets)
    gridSize = np.array(gridSize, dtype=np.int32)
    BeamletMatrix = []

    for index in range(nbOfBeamlets):
        scenarioShiftCorrected_voxel = np.round(correctShift(scenarioShift_voxel, beamletAngles_rad[index]),3)
        # scenarioShiftCorrected_voxel = scenarioShift_voxel
        beamlet = sparseBeamlets[:, index]
        shiftTrunctated = CCCdoseEngineIO.convertTo1DcoordFortran(np.trunc(scenarioShiftCorrected_voxel), gridSize) 
        nonZeroIndexes = beamlet.nonzero()
        if len(nonZeroIndexes[0]) == 0:
            BeamletMatrix.append(beamlet)  
            continue
        nonZeroValues = np.array(beamlet[nonZeroIndexes], dtype= np.float32)
        nonZeroIndexes = np.array(nonZeroIndexes[0] + shiftTrunctated, dtype= np.int32)
        NumberOfElements = np.int32(len(nonZeroValues[0]))

        nonZeroValues_gpu = cuda.mem_alloc(nonZeroValues.nbytes)
        nonZeroIndexes_gpu = cuda.mem_alloc(nonZeroIndexes.nbytes)
        gridSize_gpu = cuda.mem_alloc(gridSize.nbytes)
        cuda.memcpy_htod(nonZeroValues_gpu, nonZeroValues)
        cuda.memcpy_htod(nonZeroIndexes_gpu, nonZeroIndexes)
        cuda.memcpy_htod(gridSize_gpu, gridSize)
        # gpuarray.zeros()
        nonZeroValuesShiftedExtendedArray = np.array([])
        nonZeroIndexesShiftedExtendedArray  = np.array([])
        scenarioShiftCorrected_voxel -= np.trunc(scenarioShiftCorrected_voxel)
        if np.sum(scenarioShiftCorrected_voxel)!=0:
            for i, shift in enumerate(scenarioShiftCorrected_voxel): #### Think on puting this into the cuda code
                if shift == 0:
                    continue
                weight = np.abs(shift) / np.sum(np.abs(scenarioShiftCorrected_voxel))
                directionalShiftCorrected_voxel = np.zeros(3)
                directionalShiftCorrected_voxel[i] = shift
                magnitude = CCCdoseEngineIO.convertTo1DcoordFortran(np.trunc(directionalShiftCorrected_voxel), gridSize) 
                direction = CCCdoseEngineIO.convertTo1DcoordFortran(np.sign(directionalShiftCorrected_voxel), gridSize) 
                directionNeg = CCCdoseEngineIO.convertTo1DcoordFortran(-1 * np.sign(directionalShiftCorrected_voxel), gridSize)
                shiftValue = np.float32((directionalShiftCorrected_voxel[np.nonzero(directionalShiftCorrected_voxel)[0][0]])%1)

                nonZeroValuesShifted = np.zeros(NumberOfElements * 2).astype(np.float32)
                nonZeroValuesShifted_gpu = cuda.mem_alloc(nonZeroValuesShifted.nbytes)
                cuda.memcpy_htod(nonZeroValuesShifted_gpu, nonZeroValuesShifted)

                nonZeroIndexesShifted = np.zeros(NumberOfElements * 2).astype(np.int32)
                nonZeroIndexesShifted_gpu = cuda.mem_alloc(nonZeroIndexesShifted.nbytes)
                cuda.memcpy_htod(nonZeroIndexesShifted_gpu, nonZeroIndexesShifted)

                Shift_voxel = np.array([magnitude, directionNeg, direction],dtype=np.int32)
                scenarioShiftVoxel_gpu = cuda.mem_alloc(Shift_voxel.nbytes)
                cuda.memcpy_htod(scenarioShiftVoxel_gpu, Shift_voxel)

                shiftSparse = mod.get_function("shiftSparse")
                shiftSparse(nonZeroValues_gpu,nonZeroValuesShifted_gpu,nonZeroIndexes_gpu,nonZeroIndexesShifted_gpu, scenarioShiftVoxel_gpu, gridSize_gpu, shiftValue, NumberOfElements, grid=(int(NumberOfElements/1024)+1,1),block=(1024,1,1))

                cuda.memcpy_dtoh(nonZeroValuesShifted,nonZeroValuesShifted_gpu)
                cuda.memcpy_dtoh(nonZeroIndexesShifted,nonZeroIndexesShifted_gpu)

                nonZeroValuesShiftedExtendedArray = np.append(nonZeroValuesShiftedExtendedArray,nonZeroValuesShifted * weight)
                nonZeroIndexesShiftedExtendedArray = np.append(nonZeroIndexesShiftedExtendedArray,nonZeroIndexesShifted)
                spp = sp.csc_matrix((nonZeroValuesShiftedExtendedArray, (nonZeroIndexesShiftedExtendedArray, np.zeros(nonZeroIndexesShiftedExtendedArray.size))), shape=(nbOfVoxelInImage, 1),dtype=np.float32)
                # print('Here')
                # printaux(spp.nonzero()[0], spp[spp.nonzero()][0])
            beamlet = sp.csc_matrix((nonZeroValuesShiftedExtendedArray, (nonZeroIndexesShiftedExtendedArray, np.zeros(nonZeroIndexesShiftedExtendedArray.size))), shape=(nbOfVoxelInImage, 1),dtype=np.float32)    
        else:
            beamlet = sp.csc_matrix((nonZeroValues[0], (nonZeroIndexes, np.zeros(nonZeroIndexes.size))), shape=(nbOfVoxelInImage, 1),dtype=np.float32)      
        BeamletMatrix.append(beamlet)  
    return sp.hstack(BeamletMatrix)

def printaux(indexes, values):
    values = np.ravel(values)
    for index, value in zip(indexes, values):
        print(f"{index:<10} {value:<10}")

def find_change_indices(arr):
    arr = np.array(arr)
    # Create a boolean array where changes occur
    changes = arr[1:] != arr[:-1]
    # Use np.where to find the indices where changes occur, add 1 because we compare with the previous element
    change_indices = np.where(changes)[0] + 1
    return change_indices.tolist()

def shiftBeamlets1(sparseBeamlets, gridSize,  scenarioShift_voxel, beamletAngles_rad):
    lib = ctypes.cdll.LoadLibrary('/home/luciano/Codes/newOpenTPS/opentps/opentps_core/opentps/core/processing/doseCalculation/photons/shiftBeamlets.so')
    
    # Define the argument types and return types for the C++ function
    lib.shiftBeamlets.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
        ctypes.c_float,
        ctypes.c_int,
        ctypes.c_int
    ]
    numThreads = os.cpu_count()
    scenarioShift_voxel[2]*=-1 ### To have the setup error in LPS. Check because some signs problem
    scenarioShift_voxel[1]*=-1 ### To have the setup error in LPS. Check because some signs problem
    nbOfBeamlets = sparseBeamlets.shape[1]
    nbOfVoxelInImage = sparseBeamlets.shape[0]
    assert(len(beamletAngles_rad), nbOfBeamlets)
    gridSize = np.array(gridSize, dtype=np.int32)
    BeamletMatrix = []
    
    nonZeroIndexes = sparseBeamlets.nonzero()
    nonZeroValues = np.array(sparseBeamlets[nonZeroIndexes], dtype= np.float32)[0]
    nonZeroIndexes_beamlet = np.array(nonZeroIndexes[0], dtype= np.int32)
    indexes_beamlet = np.array(nonZeroIndexes[1], dtype= np.int32)
    arg = np.argsort(indexes_beamlet)
    
    indexes_beamlet = indexes_beamlet[arg]
    nonZeroIndexes_beamlet = nonZeroIndexes_beamlet[arg]
    nonZeroValues = nonZeroValues[arg]
    indexes = [0] + find_change_indices(indexes_beamlet)
    
    # lib.shiftBeamlets(nonZeroValues, nonZeroValuesShifted, nonZeroIndexes, nonZeroIndexesShifted, Shift_voxel, shiftValue, NumberOfElements, numThreads)   
        
    return sp.hstack(BeamletMatrix)


def shiftBeamletscpp(sparseBeamlets, gridSize,  scenarioShift_voxel, beamletAngles_rad):
    lib = ctypes.cdll.LoadLibrary('/home/luciano/Codes/newOpenTPS/opentps/opentps_core/opentps/core/processing/doseCalculation/photons/shiftBeamlets.so')
    
    # Define the argument types and return types for the C++ function
    lib.shiftBeamlets.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        ctypes.c_int,
        ctypes.c_int
    ]
    numThreads = psutil.cpu_count()
    scenarioShift_voxel[2]*=-1 ### To have the setup error in LPS. Check because some signs problem
    scenarioShift_voxel[1]*=-1 ### To have the setup error in LPS. Check because some signs problem
    nbOfBeamlets = sparseBeamlets.shape[1]
    nbOfVoxelInImage = sparseBeamlets.shape[0]
    assert(len(beamletAngles_rad), nbOfBeamlets)
    gridSize = np.array(gridSize, dtype=np.int32)
    BeamletMatrix = []
    
    length = sparseBeamlets.shape[0]
    nonZeroIndexes = sparseBeamlets.nonzero()
    nonZeroValues = np.array(sparseBeamlets[nonZeroIndexes], dtype= np.float32)[0]
    nonZeroIndexes_beamlet = np.array(nonZeroIndexes[0], dtype= np.int32)
    indexes_beamlet = np.array(nonZeroIndexes[1], dtype= np.int32)
    arg = np.argsort(indexes_beamlet)
    
    indexes_beamlet = indexes_beamlet[arg]
    nonZeroIndexes_beamlet = nonZeroIndexes_beamlet[arg]
    nonZeroValues = nonZeroValues[arg]
    indexesChangeBeamlet = [0] + find_change_indices(indexes_beamlet) + [len(indexes_beamlet)]
    NumberOfElements = len(nonZeroValues)
    nonZeroValues = np.array(nonZeroValues, dtype=np.float32)
    nonZeroIndexes_beamlet = np.array(nonZeroIndexes_beamlet, dtype=np.int32)
    indexesChangeBeamlet = np.array(indexesChangeBeamlet, dtype=np.int32)
    scenarioShift_voxel = np.array(scenarioShift_voxel, dtype=np.float32)
    
    nonZeroValuesShifted = np.zeros(NumberOfElements * 2 * 3).astype(np.float32)
    nonZeroIndexesShifted = np.zeros(NumberOfElements * 2 * 3).astype(np.int32)
    beamletAngles_rad = np.array(beamletAngles_rad, dtype=np.float32)
    nOfBeamlets = len(beamletAngles_rad)
    lib.shiftBeamlets(nonZeroValues, nonZeroValuesShifted, nonZeroIndexes_beamlet, nonZeroIndexesShifted, indexesChangeBeamlet, scenarioShift_voxel, gridSize, beamletAngles_rad, nOfBeamlets, numThreads)
    
    for i in range(nOfBeamlets):
        start = indexesChangeBeamlet[i]
        end = indexesChangeBeamlet[i+1]
        indexes = nonZeroIndexesShifted[start*2*3:end*2*3]
        values = nonZeroValuesShifted[start*2*3:end*2*3]
        indexes = indexes[indexes.nonzero()]
        values = values[values.nonzero()]
        BeamletMatrix.append(sp.csc_matrix((values, (indexes, np.zeros(len(indexes)))), shape=(length,1), dtype=np.float32))
        
    return sp.hstack(BeamletMatrix, format='csc')

def convolveSparseMatrix(sparse, sigma_voxels, image_size):
    len = sparse.shape[0]
    matrix = np.reshape(sparse.A, image_size ,order='F')
    convolved_matrix = ndimage.gaussian_filter(matrix.astype(float), sigma = sigma_voxels, order=0, truncate=2)
    matrix = np.reshape(convolved_matrix, (len, 1), order='F')
    return sp.csc_matrix(matrix)


def dnorm(x, mu, sd):
    return 1 / (np.sqrt(2 * np.pi) * sd) * np.e ** (-np.power((x - mu) / sd, 2) / 2)

def gaussian_kernel_3d(size, sigma=1): ### Gaussian Kernel used to smooth the result of the sigma smooth
    if sigma==0:
        sigma = 1e-6
    kernel_1D = np.linspace(-(size // 2), size // 2, size)
    for i in range(size):
        kernel_1D[i] = dnorm(kernel_1D[i], 0, sigma)
    kernel_2D = np.outer(kernel_1D.T, kernel_1D.T)
    kernel_3D = np.einsum("ij,k->ijk", kernel_2D,kernel_1D)
    kernel_3D = kernel_3D / np.sum(kernel_3D)
    return kernel_3D


def adjustDoseToScenario(scenario, nominal, imageSpacing, plan: PhotonPlan): #### it might not fit here
    if scenario.sse is not None:
        shiftVoxels = np.array(scenario.sse) / np.array(imageSpacing)
        cumulativeNumberBeamlets = 0
        weights = nominal.sb.beamletWeights
        doseGridSize = nominal.sb.doseGridSize
        dose = np.zeros(doseGridSize)
        sizeImage = nominal.sb._sparseBeamlets.shape[0]
        nofBeamlets = nominal.sb._sparseBeamlets.shape[1]
        assert(nofBeamlets==len(plan.beamlets), f"The number of beamlets in the dose influece matrix is {nofBeamlets} but the number of beamlets in the treatment plan is {len(plan.beamlets)}")
        for segment in plan.beamSegments:
            beamletsSegment = nominal.sb._sparseBeamlets[:, cumulativeNumberBeamlets: cumulativeNumberBeamlets + len(segment)]
            weightsSegment = weights[cumulativeNumberBeamlets: cumulativeNumberBeamlets + len(segment)]
            result = csc_matrix.dot(beamletsSegment, weightsSegment).reshape(sizeImage,1)
            result = np.reshape(result, doseGridSize, order='F')
            result = np.flip(result, 0)
            result = np.flip(result, 1)
            shiftVoxelsCorrected = np.round(correctShift(shiftVoxels, segment.gantryAngle_degree / 180 * np.pi),3)
            dose +=  shift(result, shiftVoxelsCorrected, mode='constant', cval=0, order=1)
            cumulativeNumberBeamlets+=len(segment)

        dose = DoseImage(imageArray=dose, origin=nominal.sb.doseOrigin, spacing=nominal.sb.doseSpacing,
                              angles=(1, 0, 0, 0, 1, 0, 0, 0, 1))
    else:
        dose = nominal.sb.toDoseImage()

    doseArray = dose.imageArray
    if scenario.sre != None:
        doseArray = gaussian_filter(doseArray.astype(float), sigma = scenario.sre, order=0, truncate=2)
    else:
        return dose
    dose.imageArray = doseArray

    return dose