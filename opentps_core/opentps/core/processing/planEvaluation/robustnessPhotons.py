from __future__ import annotations

import logging
from enum import Enum
from typing import Union

import numpy as np
import pickle
import os

from opentps.core.data._sparseBeamlets import SparseBeamlets
import copy
from opentps.core.data._dvhBand import DVHBand
from opentps.core.processing.imageProcessing import resampler3D
from opentps.core.data.images import ProbabilityMap, cleanArrayWithThreshold

from typing import TYPE_CHECKING

# if TYPE_CHECKING:
from opentps.core.data import ROIContour
from opentps.core.data.images import ROIMask, DoseImage

logger = logging.getLogger(__name__)
from opentps.core.data._dvh import DVH


class RobustScenario:
    def __init__(self, sb:SparseBeamlets = None, sse = None, sre = None, pdf = 0, dilation_mm = {}):
        self.sse = sse      # Setup Systematic Error
        self.sre = sre      # Setup Random Error
        self.sb = sb
        self.dilation_mm = dilation_mm
        self.pdf = pdf
        self.dose = []
        self.dvh = []
        self.targetD95 = 0
        self.targetD5 = 0
        self.targetMSE = 0
        self.selected = 0

    def printInfo(self):
        logger.info('Setup Systematic Error:{} mm'.format(self.sse))
        logger.info('Setup Random Error:{} mm'.format(self.sre))
        logger.info("Target_D95 = " + str(self.targetD95))
        logger.info("Target_D5 = " + str(self.targetD5))
        logger.info("Target_MSE = " + str(self.targetMSE)+"\n")
        # logger.info(" ")

    def save(self, file_path):
        with open(file_path, 'wb') as fid:
            pickle.dump(self.__dict__, fid)
    
    def copy(self):
        return copy.deepcopy(self)

    def load(self, file_path):
        with open(file_path, 'rb') as fid:
            tmp = pickle.load(fid)

        self.__dict__.update(tmp)

    def unload(self):
        if self.sb != None:
            self.sb.unload()
    
    def toSparseMatrix(self):
        return self.sb.toSparseMatrix()
    
    def setUnitaryBeamlets(self,beamletMatrix):
        self.sb.setUnitaryBeamlets(beamletMatrix)

    def __str__(self):
        str = 'Setup Systematic Error:{} mm\n'.format(self.sse) + 'Setup Random Error:{} mm\n'.format(self.sre)
        for name, dilation in self.dilation_mm.items():
            str += name + f' dilated {round(dilation,2)} mm \n'
        return str
    
    def toTxt(self, path):
        f = open(path, "w")
        f.write(str(self))
        f.close()

class Robustness:

    class Strategies(Enum): #### Use these classes
        DEFAULT = "DISABLED"
        DISABLED = "DISABLED"
        ERRORSPACE_REGULAR = "ERRORSPACE_REGULAR"
        ERRORSPACE_STAT = "ERRORSPACE_STAT"
        DOSIMETRIC = "DOSIMETRIC"
    
    class Modality(Enum):
        MINMAX = "MINMAX"
        PRO = "PRO"

    def __init__(self):
        self.selectionStrategy = self.Strategies.DEFAULT
        self.setupSystematicError = None  # mm
        self.numberOfSigmas = 2.5 # Number of sigmas in the normal distribution
        self.setupRandomError = None  # mm
        self.target = []
        self.targetPrescription = 60  # Gy
        self.nominal = RobustScenario()
        # self.numScenarios = 0
        self.scenarios:list[RobustScenario] = []
        self.dvhBands = []
        self.doseDistributionType = ""
        self.doseDistribution = []
        self.sseNumberOfSamples = 1
        self.modality = self.Modality.MINMAX
        

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index == 0:
            self.index += 1
            return self.nominal
        if self.index < len(self.scenarios)+1:
            result = self.scenarios[self.index - 1]
            self.index += 1
            return result
        else:
            raise StopIteration
        
    def setNominal(self, sb:SparseBeamlets):
        self.nominal.sse = [0,0,0]
        self.nominal.sb = sb
        self.nominal.pdf = self.pdf(self.nominal.sse)

    def pdf(self, vector):
        x,y,z = vector
        var = np.array(self.setupSystematicError)**2
        return np.exp(-1/2*(x**2 / var[0] + y**2 /var[1] + z**2 /var[2])) / np.sqrt((2*np.pi)**3 * np.prod(var))

    def addScenario(self, sb:SparseBeamlets, sse , sre):
        scenario = RobustScenario(sb, sse , sre, self.pdf(sse))
        self.scenarios.append(scenario)
    
    def generateRobustScenarios4Planning(self):
        if self.setupSystematicError!=None:
            for index, sse in enumerate(self.setupSystematicError):
                for scale in range(1,self.sseNumberOfSamples+1):
                    sseScaled = sse / self.sseNumberOfSamples * scale
                    for sign in [-1,1]:
                        array = np.zeros(3)
                        array[index] = sseScaled * sign * self.numberOfSigmas
                        scenario = RobustScenario(sse = array, pdf = self.pdf(array))
                        self.scenarios.append(scenario)



    def sampleScenario(self):
        sse_sampled = list(np.random.normal([0]*len(self.setupSystematicError),self.setupSystematicError)) if self.setupSystematicError != None else None
        sre_sampled = np.random.uniform([self.setupRandomError,0]) if self.setupRandomError != None else None
        return RobustScenario(sb = None, sse = sse_sampled, sre = sre_sampled)
    
    def unload(self):
        self.nominal.unload()
        for scenario in self.scenarios:
            scenario.unload()

    def setBeamWeights(self, spotMUs):
        self.nominal.sb.beamletWeights = spotMUs
        for scenario in self.scenarios:
            scenario.sb.beamletWeights = spotMUs
    
    @property
    def numScenarios(self): ### Should be only RO scenarios or nominal as well?
        return len(self.scenarios) 

    @property
    def pdfScenarios(self): 
        return np.array([scenario.pdf for scenario in self])

    def load(self, path):
        with open(path, 'rb') as file:
            tmp = pickle.load(file)
        self.__dict__.update(tmp)
    
    def save(self, path):
        with open(path, 'wb') as file:
            pickle.dump(self, file)


class RobustEvaluation(Robustness):
    
    def __init__(self, robustness = None, contours = None):
        if robustness != None and contours != None:
            self.copy_attributes(robustness)
            contours = self.convertToNominal(contours)
            if self.nominal.sb != None:
                self.initializeDVHs(self.nominal, self.nominal.sb.toDoseImage(), contours=contours)
            if len(self.scenarios)>0:
                for scenario in self.scenarios:
                    self.initializeDVHs(scenario, scenario.sb.toDoseImage(), contours=contours)
            self.unload()
        else:
            super().__init__()

    def copy_attributes(self, robustness):
        for attr_name, attr_value in vars(robustness).items():
            setattr(self, attr_name, copy.deepcopy(attr_value))

    def convertToNominal(self, structures):
        samples = []
        for structure in structures:
            if isinstance(structure, ProbabilityMap):
                samples.append(structure.getCTV())
            else:
                samples.append(structure)
            assert(isinstance(samples[-1], ROIMask))
        return samples

    def initializeDVHs(self, scenario, dose: DoseImage, contours: Union[ROIContour, ROIMask]):
        scenario.dose = dose
        # Need to set patient to None for memory, est-ce que ca va poser probleme ?
        scenario.dose.patient = None
        scenario.dvh.clear()
        for contour in contours:
            contour.patient = None
            myDVH = DVH(contour, scenario.dose)
            scenario.dvh.append(myDVH)
        scenario.dose.imageArray = scenario.dose.imageArray.astype(np.float32)  # can be reduced to float16 because all metrics are already computed and it's only used for display

    def setNominal(self, dose: DoseImage, contours: Union[ROIContour, ROIMask]):
        self.initializeDVHs(self.nominal, dose, contours=contours)

    def addScenario(self, dose: DoseImage, contours: Union[ROIContour, ROIMask], scenario = None):
        if scenario == None: 
            scenario = RobustScenario()
        self.initializeDVHs(scenario, dose, contours=contours)
        scenario.unload()
        self.scenarios.append(scenario)

    def setTarget(self, targetContour, targetPrescription):
        if not(self.nominal.dose.hasSameGrid(targetContour)):
            resampler3D.resampleImage3DOnImage3D(targetContour,self.nominal.dose, inPlace=True, fillValue=0.)
        self.target = targetContour
        self.target.imageArray = self.target.imageArray > 0
        self.targetPrescription = targetPrescription
        for dvh in self.nominal.dvh:
            if dvh._roiName == self.target.name:
                self.nominal.targetD95 = dvh.D95
                self.nominal.targetD5 = dvh.D5
                self.nominal.targetMSE = self.computeTargetMSE(self.nominal.dose.imageArray)
                break

        for scenario in self.scenarios:
            for dvh in scenario.dvh:
                if dvh._roiName == self.target.name:
                    scenario.targetD95 = dvh.D95
                    scenario.targetD5 = dvh.D5
                    scenario.targetMSE = self.computeTargetMSE(scenario.dose.imageArray)
                    break

    def recomputeDVH(self, contours):
        from opentps.core.data._dvh import DVH
        self.nominal.dvh.clear()
        for contour in contours:
            myDVH = DVH(contour, self.nominal.dose)
            self.nominal.dvh.append(myDVH)

        for scenario in self.scenarios:
            scenario.dvh.clear()
            for contour in contours:
                myDVH = DVH(contour, scenario.dose)
                scenario.dvh.append(myDVH)

    def computeTargetMSE(self, dose):
        dose_vector = dose[self.target.imageArray]
        error = dose_vector - self.targetPrescription
        mse = np.mean(np.square(error))
        return mse

    def analyzeErrorSpace(self, metric, targetContour, targetPrescription):
        if (
                self.target == [] or self.target.name != targetContour.name or self.targetPrescription != targetPrescription):
            self.setTarget(targetContour, targetPrescription)

        # sort scenarios from worst to best according to selected metric
        if metric == "D95":
            self.scenarios.sort(key=(lambda scenario: scenario.targetD95))
        elif metric == "MSE":
            self.scenarios.sort(key=(lambda scenario: scenario.targetMSE))

        # initialize dose distribution
        if self.doseDistributionType == "Nominal":
            self.doseDistribution = self.nominal.dose.copy()
        else:
            self.doseDistribution = self.scenarios[0].dose.copy()  # Worst scenario

        # initialize dvh-band structure
        allDVH = []
        allDmean = []
        for dvh in self.scenarios[0].dvh:
            allDVH.append(np.array([]).reshape((len(dvh._volume), 0)))
            allDmean.append([])

        # generate DVH-band
        for s in range(self.numScenarios):
            self.scenarios[s].selected = 1
            if self.doseDistributionType == "Voxel wise minimum":
                self.doseDistribution.imageArray = np.minimum(self.doseDistribution.imageArray, self.scenarios[s].dose.imageArray)
            elif self.doseDistributionType == "Voxel wise maximum":
                self.doseDistribution.imageArray = np.maximum(self.doseDistribution.imageArray, self.scenarios[s].dose.imageArray)
            for c in range(len(self.scenarios[s].dvh)):
                allDVH[c] = np.hstack((allDVH[c], np.expand_dims(self.scenarios[s].dvh[c]._volume, axis=1)))
                allDmean[c].append(self.scenarios[s].dvh[c].Dmean)

        self.dvhBands.clear()
        for c in range(len(self.scenarios[0].dvh)):
            dvh = self.scenarios[0].dvh[c]
            dvhBand = DVHBand()
            dvhBand._roiName = dvh._roiName
            dvhBand._dose = dvh._dose
            dvhBand._volumeLow = np.amin(allDVH[c], axis=1)
            dvhBand._volumeHigh = np.amax(allDVH[c], axis=1)
            dvhBand._nominalDVH = self.nominal.dvh[c]
            dvhBand.computeMetrics()
            dvhBand._Dmean = [min(allDmean[c]), max(allDmean[c])]
            self.dvhBands.append(dvhBand)

    def analyzeDosimetricSpace(self, metric, CI, targetContour, targetPrescription):
        if (
                self.target == [] or self.target.name != targetContour.name or self.targetPrescription != targetPrescription):
            self.setTarget(targetContour, targetPrescription)

        if metric == "D95":
            self.scenarios.sort(key=(lambda scenario: scenario.targetD95))
        elif metric == "MSE":
            self.scenarios.sort(key=(lambda scenario: scenario.targetMSE))

        start = round(self.numScenarios * (100 - CI) / 100)
        if start == self.numScenarios: start -= 1

        # initialize dose distribution
        if self.doseDistributionType == "Nominal":
            self.doseDistribution = self.nominal.dose.copy()
        else:
            self.doseDistribution = self.scenarios[start].dose.copy()  # Worst scenario

        # initialize dvh-band structure
        selectedDVH = []
        selectedDmean = []
        for dvh in self.scenarios[0].dvh:
            selectedDVH.append(np.array([]).reshape((len(dvh.volume), 0)))
            selectedDmean.append([])

        # select scenarios
        for s in range(self.numScenarios):
            if s < start:
                self.scenarios[s].selected = 0
            else:
                self.scenarios[s].selected = 1
                if self.doseDistributionType == "Voxel wise minimum":
                    self.doseDistribution.imageArray = np.minimum(self.doseDistribution.imageArray, self.scenarios[s].dose.imageArray)
                elif self.doseDistributionType == "Voxel wise maximum":
                    self.doseDistribution.imageArray = np.maximum(self.doseDistribution.imageArray, self.scenarios[s].dose.imageArray)
                for c in range(len(self.scenarios[s].dvh)):
                    selectedDVH[c] = np.hstack(
                        (selectedDVH[c], np.expand_dims(self.scenarios[s].dvh[c].volume, axis=1)))
                    selectedDmean[c].append(self.scenarios[s].dvh[c].Dmean)

        # compute DVH-band envelopes
        self.dvhBands.clear()
        for c in range(len(self.scenarios[s].dvh)):
            dvh = self.scenarios[0].dvh[c]
            dvhBand = DVHBand()
            dvhBand._roiName = dvh._roiName
            dvhBand._dose = dvh._dose
            dvhBand._volumeLow = np.amin(selectedDVH[c], axis=1)
            dvhBand._volumeHigh = np.amax(selectedDVH[c], axis=1)
            dvhBand._nominalDVH = self.nominal.dvh[c]
            dvhBand.computeMetrics()
            dvhBand._Dmean = [min(selectedDmean[c]), max(selectedDmean[c])]
            self.dvhBands.append(dvhBand)

    def printInfo(self):
        logger.info("Nominal scenario:")
        self.nominal.printInfo()

        for i in range(len(self.scenarios)):
            logger.info("Scenario " + str(i + 1))
            self.scenarios[i].printInfo()

    def save(self, folder_path):
        try:
            file_path = os.path.join(folder_path, "RobustnessEvaluation" + ".pkl")
            with open(file_path, "wb") as file:
                pickle.dump(self, file)
        except:
            print('could not save the RO error happened')

    def load(self, folder_path):
        file_path = os.path.join(folder_path, "RobustnessEvaluation" + ".pkl")
        with open(file_path, "rb") as fid:
            tmp = pickle.load(fid)
        self.__dict__.update(tmp.__dict__)


