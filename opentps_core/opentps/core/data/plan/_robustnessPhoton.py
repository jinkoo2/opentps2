__all__ = ['RobustnessPhoton','RobustScenario']

import copy
import itertools
import pickle
import logging
import numpy as np
from opentps.core.data._sparseBeamlets import SparseBeamlets
from opentps.core.data.plan._robustness import Robustness

logger = logging.getLogger(__name__)


class RobustScenario:
    def __init__(self, sb:SparseBeamlets = None, sse = None, sre = None, dilation_mm = {}):
        self.sse = sse      # Setup Systematic Error
        self.sre = sre      # Setup Random Error
        self.sb = sb
        self.dilation_mm = dilation_mm
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

class RobustnessPhoton(Robustness):
    """
    This class creates an object that stores the robustness parameters of a photon plan and generates robust scenarios through sampling(optimization).

    Attributes
    ----------
    rangeSystematicError : float (default = 1.6) (%)
        The range systematic error in %.
    scenarios : list
        The list of scenarios.
    """

    def __init__(self):
        self.numberOfSigmas = 2.5
        self.sseNumberOfSamples = 1
        self.scenarios:list[RobustScenario] = []
        self.nominal = RobustScenario()
        
        
        super().__init__()

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

    def addScenario(self, sb:SparseBeamlets, sse , sre):
        scenario = RobustScenario(sb, sse , sre)
        self.scenarios.append(scenario)
    
    def generateRobustScenarios(self):
        if self.setupSystematicError not in [None, 0, [0,0,0]]:
            if self.selectionStrategy == self.selectionStrategy.RANDOM :
                self.generateRandomScenarios()
            elif self.selectionStrategy == self.selectionStrategy.REDUCED_SET :
                self.generateReducedErrorSpacecenarios()
            elif self.selectionStrategy == self.selectionStrategy.ALL :
                self.generateAllErrorSpaceScenarios()
            elif self.selectionStrategy == self.selectionStrategy.DISABLED :
                self.scenarios.append(RobustScenario(sse = np.array([self.setupSystematicError[0]* self.numberOfSigmas, 
                                                                         self.setupSystematicError[1]* self.numberOfSigmas, 
                                                                         self.setupSystematicError[2]* self.numberOfSigmas])))
        else :
            raise Exception("No evaluation strategy selected")

    def generateReducedErrorSpacecenarios(self):  # From [a, b, c] to 6 scenarios [+-a, +-b, +-c]
        for index, sse in enumerate(self.setupSystematicError):
            for sign in [-1,1]:
                array = np.zeros(3)
                array[index] = sse * sign * self.numberOfSigmas
                scenario = RobustScenario(sse = array, sre = self.setupRandomError)
                self.scenarios.append(scenario)

    def generateAllErrorSpaceScenarios(self):
        # Point coordinates on hypersphere with two zero axes
        R = self.setupSystematicError[0] * self.numberOfSigmas
        for sign in [-1, 1]:
            self.scenarios.append(RobustScenario(sse = np.round(np.array([sign * R, 0, 0]), 2), sre = self.setupRandomError))
            self.scenarios.append(RobustScenario(sse = np.round(np.array([0, sign * R, 0]), 2), sre = self.setupRandomError))
            self.scenarios.append(RobustScenario(sse = np.round(np.array([0, 0, sign * R]), 2), sre = self.setupRandomError))

        # Coordinates of point on hypersphere with zero axis
        sqrt2 = R / np.sqrt(2)
        for sign1, sign2 in itertools.product([-1, 1], repeat=2):
            self.scenarios.append(RobustScenario(sse = np.round(np.array([sign1 * sqrt2, sign2 * sqrt2, 0]), 2), sre = self.setupRandomError))
            self.scenarios.append(RobustScenario(sse = np.round(np.array([sign1 * sqrt2, 0, sign2 * sqrt2]), 2), sre = self.setupRandomError))
            self.scenarios.append(RobustScenario(sse = np.round(np.array([0, sign1 * sqrt2, sign2 * sqrt2]), 2), sre = self.setupRandomError))

        # Coordinates of point on hypersphere without any zero axis (diagonals)
        sqrt3 = R / np.sqrt(3)
        for signs in itertools.product([-1, 1], repeat=3):
            self.scenarios.append(RobustScenario(sse = np.round(np.array([signs[0] * sqrt3, signs[1] * sqrt3, signs[2] * sqrt3]), 2), sre = self.setupRandomError))


    def generateRandomScenarios(self):
        # Sample in gaussian
        setupErrorSpace = self.setupSystematicError
        for _ in range(self.NumScenarios):
            SampleSetupError = [np.random.normal(0, sigma) for sigma in setupErrorSpace]
            SampleSetupError = np.round(SampleSetupError, 2)
            scenario = RobustScenario(sse = SampleSetupError, sre = self.setupRandomError)
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

    def load(self, path):
        with open(path, 'rb') as file:
            tmp = pickle.load(file)
        self.__dict__.update(tmp)
    
    def save(self, path):
        with open(path, 'wb') as file:
            pickle.dump(self, file)