__all__ = ['ProbabilityMap', 'cleanArrayWithThreshold']

import pydicom
import copy

from opentps.core.data.images._roiMask import ROIMask
from opentps.core.io.sitkIO import exportImageSitk
from opentps.core.processing.segmentation.segmentationModel import ProbabilisticModel
import numpy as np
import scipy.stats as stats
from scipy.ndimage import morphology
import os
from scipy import ndimage

# self.imageArray has the binary contour that correspond to the whole extention of the probability map
class ProbabilityMap(ROIMask): 
    def __init__(self, image, probabilisticModel: ProbabilisticModel = None, probabilityThreshold = 0.01):
        extensionName = ''
        self.probabilityThreshold = probabilityThreshold
        self.probabilisticModel = probabilisticModel
        if probabilisticModel is not None:
            self.gtv = copy.deepcopy(image)
            ROIMask.__init__(self, imageArray=image.imageArray, name=image.name, origin=image.origin, spacing=image.spacing, angles=image.angles, seriesInstanceUID=image.seriesInstanceUID, patient= image.patient)    
            self.probabilityMap = probabilisticModel(image)
        else:
            self.probabilityMap = image.imageArray
            self.probabilityMap[self.probabilityMap < self.probabilityThreshold] = 0
            ROIMask.__init__(self, imageArray = image.imageArray > self.probabilityThreshold, name=image.name, origin=image.origin, spacing=image.spacing, angles=image.angles,
                            seriesInstanceUID=image.seriesInstanceUID, patient= image.patient)
            # self.name = "Probability Map_" + image.name
        self.updateImageArray() 

    def __str__(self):
        return "Probability Map: " + self.seriesInstanceUID

    def thresholdProbMap(self,thr):
        return ROIMask(imageArray= self.probabilityMap>=thr, name= str(thr) + '_' + self.name, origin=self.origin, spacing=self.spacing, angles=self.angles, seriesInstanceUID=self.seriesInstanceUID, patient= self.patient)

    def cleanWithThreshold(self, ct, threshold):
        self.probabilityMap = cleanArrayWithThreshold(self.probabilityMap, ct, threshold) 
        self.updateImageArray()
        self.probabilityMap *= largest_component(self.imageArray)
        self.updateImageArray()

    def setProbabilityMap(self, pMap):
        self.probabilityMap = pMap
        self.updateImageArray()

    def updateImageArray(self):
        self.imageArray = self.probabilityMap > self.probabilityThreshold

    def getCTV(self,probMapThr = 0.9):
        array = self.thresholdProbMap(probMapThr).imageArray
        return ROIMask(imageArray = array, name = self.name, origin=self.origin, spacing=self.spacing, angles=self.angles,
                         seriesInstanceUID=self.seriesInstanceUID, patient= self.patient)
    
    def sample(self): ### Double check
        if self.probabilisticModel is not None:
            dilate_mm = np.sum(self.probabilisticModel.sample())
            ctvSampled = copy.deepcopy(self.gtv)
            ctvSampled.dilate(dilate_mm)
            return ctvSampled, dilate_mm
        else:
            x,y,z = self.probabilityMap.shape
            sample = (np.random.rand(x,y,z) < self.probabilityMap).astype(int)
            return ROIMask(imageArray = sample, name = self.name, origin=self.origin, spacing=self.spacing, angles=self.angles,
                            seriesInstanceUID=self.seriesInstanceUID, patient= self.patient)

    def sampleFromRefMask(self,refMask):
        sample = self.probabilisticModel.sample()
        dilate_mm = sample if not isinstance(sample, list) else sample[-1]
        ctvSampled = copy.deepcopy(refMask)
        ctvSampled.dilate(dilate_mm)
        return ctvSampled, dilate_mm    

    def saveProbMap(self,path):
        if not os.path.isdir(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        image = self.copy()
        image.imageArray = self.probabilityMap
        exportImageSitk(path, image)

    def saveProbMapMask(self,path):
        image = self.copy()
        image.imageArray = self.imageArray * 1
        exportImageSitk(path, image)  

    def dilate(self, dilate_mm):
        if dilate_mm>0:
            if self.probabilisticModel is not None:
                # gtvCopy = copy.deepcopy(self.gtv)
                self.gtv.dilate(dilate_mm) ### Not sure
                self.probabilityMap = self.probabilisticModel(self.gtv)
                self.updateImageArray() 
            else:
                NotImplementedError('No implemented sampling for PM which is not CTM')

def removeFromMask(mask, CT, minThreshold, maxThreshold):
    return mask - mask * np.logical_and(CT > minThreshold,CT < maxThreshold)

def cleanArrayWithThreshold(array, ct, threshold):
    for thr in threshold:
        array = removeFromMask(array, ct.imageArray, thr["min"], thr["max"]) 
    return array


def largest_component(binary_image):
    # Label connected components in the binary image
    labeled_image, num_features = ndimage.label(binary_image)

    # Calculate the size of each labeled component
    component_sizes = np.bincount(labeled_image.flatten())

    # Find the label of the largest component (excluding background)
    largest_component_label = np.argmax(component_sizes[1:]) + 1

    # Create a binary mask for the largest component
    largest_component_mask = (labeled_image == largest_component_label)

    # Apply the mask to the original binary image
    # result_image = binary_image * largest_component_mask

    return largest_component_mask