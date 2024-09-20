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