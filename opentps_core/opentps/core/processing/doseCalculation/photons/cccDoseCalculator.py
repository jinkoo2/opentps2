import copy
import logging
import math
import os
import platform
import shutil
import subprocess
import scipy as sp
from pathlib import Path
from typing import Optional, Sequence, Union
from matplotlib import pyplot as plt
from scipy.sparse import csc_matrix


import numpy as np
from opentps.core.data.images import DoseImage
from opentps.core.data import SparseBeamlets
from opentps.core.processing.doseCalculation.abstractDoseCalculator import AbstractDoseCalculator
from opentps.core.utils.programSettings import ProgramSettings
from opentps.core.data.CTCalibrations._abstractCTCalibration import AbstractCTCalibration
from opentps.core.data.images import CTImage
from opentps.core.data.images import Image3D
from opentps.core.data.images import ROIMask
from opentps.core.data import ROIContour
from opentps.core.data.plan._photonPlan import PhotonPlan
import opentps.core.io.CCCdoseEngineIO as CCCdoseEngineIO
from opentps.core.processing.doseCalculation.photons._utils import shiftBeamlets, adjustDoseToScenario

import time
# from opentps.core.processing.planEvaluation.robustnessPhotons import Robustness

__all__ = ['CCCDoseCalculator']


logger = logging.getLogger(__name__)


class CCCDoseCalculator(AbstractDoseCalculator):
    """
    Class for Collapse Cone Convolution dose calculation algorithm using WiscPlan Engine.
    This class is a wrapper for the Collapse Cone Convolution dose calculation algorithm.

    Attributes
    ----------
    _ctCalibration : AbstractCTCalibration
        CT calibration (Optional)
    _ct : Image3D
        the CT image of the patient (Optional)
    _plan : PhotonPlan
        Treatment plan (Optional)
    _roi : ROIMask
        ROI mask
    _simulationDirectory : str
        Simulation directory path
    _simulationFolderName : str
        Simulation folder name
    batchSize : int
        number of processes created in the cpu to calculate the dose
    _subprocess : subprocess.Popen
        Subprocess if used
    _subprocessKilled : bool
        Subprocess killed (if subprocess is used)
    overwriteOutsideROI : bool
        if true, set to air all the region in the CT outside the ROI
    """
    def __init__(self, batchSize = 1):

        self._ctCalibration: Optional[AbstractCTCalibration] = None
        self._ct: Optional[Image3D] = None
        self._plan: Optional[PhotonPlan] = None
        self._roi = None
        self._simulationDirectory = ProgramSettings().simulationFolder
        self._simulationFolderName = 'CCC_simulation'
        self.batchSize = batchSize

        self._subprocess = []
        self._subprocessKilled = True

        self.overwriteOutsideROI = None  # Previously cropCTContour but this name was confusing

        current_dir = os.path.dirname(os.path.abspath(__file__))
        self.WorkSpaceDir = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir, os.pardir, os.pardir, os.pardir))
        self.ROFolder = ''
    
    @property
    def _CCCSimuDir(self):
        folder = os.path.join(self._simulationDirectory, self._simulationFolderName, self.ROFolder)
        self._createFolderIfNotExists(folder)
        return folder

    @property
    def outputDir(self):
        folder = os.path.join(self._CCCSimuDir, 'Outputs')
        self._createFolderIfNotExists(folder)
        return folder

    @property
    def _executableDir(self):
        dir = os.path.join(self._CCCSimuDir, 'execFiles')
        self._createFolderIfNotExists(dir)
        return dir

    @property
    def _ctName(self):
        return 'Geometry'
    
    @property
    def _ctDirName(self):
        ctDir = os.path.join(self._CCCSimuDir, self._ctName)
        self._createFolderIfNotExists(ctDir)
        return ctDir
    
    @property
    def _beamDirectory(self):
        dir = os.path.join(self._CCCSimuDir, 'BeamSpecs')
        self._createFolderIfNotExists(dir)
        return dir
    
    @property
    def _CCCexecutablePath(self):
        return os.path.join(self.WorkSpaceDir,'opentps','core','processing','doseCalculation','photons','CCC_DoseEngine', 'CCC_DoseEngine')
    
    @property
    def ctCalibration(self) -> Optional[AbstractCTCalibration]:
        return self._ctCalibration

    @ctCalibration.setter
    def ctCalibration(self, ctCalibration: AbstractCTCalibration):
        self._ctCalibration = ctCalibration

    def createKernelFilePath(self):
        kernelsDir = os.path.join(self.WorkSpaceDir,'opentps','core','processing','doseCalculation','photons','Kernels_differentFluence')
        f = open(os.path.join(self._CCCSimuDir, 'kernelPaths.txt'),'w')
        for fileName in os.listdir(kernelsDir):
            split = fileName.split('.')
            if split[-1] == 'txt':
                f.write(split[0]+'\n')
            else:
                f.write('kernel_'+split[0]+'\n')
            f.write(os.path.join(kernelsDir, fileName)+'\n')                

    @property
    def _kernelsFilePath(self):
        kernelFilePath = os.path.join(self._CCCSimuDir, 'kernelPaths.txt')
        # if not os.path.isfile(kernelFilePath):
        self.createKernelFilePath()
        return kernelFilePath

    def createGeometryFilePath(self):
        f = open(os.path.join(self._CCCSimuDir, 'geometryFilePath.txt'),'w')
        f.write('geometry_header\n'+os.path.join(self._ctDirName, 'CT_HeaderFile.txt\n'))
        f.write('geometry_density\n'+os.path.join(self._ctDirName, 'CT.bin\n'))
  

    @property
    def _geometryFilePath(self):
        geometryFilePath = os.path.join(self._CCCSimuDir, 'geometryFilePath.txt')
        if not os.path.isfile(geometryFilePath):
            self.createGeometryFilePath()
        return geometryFilePath

    
    def writeExecuteCCCfile(self):
        for batch in range(self.batchSize):
            f = open(os.path.join(self._executableDir, 'CCC_simulation_batch{}'.format(batch)),'w')
            f.write('{executablePath} {kernelFilePath} {geometryFilePath} {beamPath} {outputPath}'.format(executablePath = self._CCCexecutablePath, kernelFilePath = self._kernelsFilePath, geometryFilePath = self._geometryFilePath, beamPath = os.path.join(self._beamDirectory,'pencilBeamSpecs_batch{}.txt'.format(batch)), outputPath = os.path.join(self.outputDir,'sparseBeamletMatrix_batch{}.bin'.format(batch))))
            f.close()


    def _writeFilesToSimuDir(self):
        CCCdoseEngineIO.writeCT(self._ct, self._ctDirName, self._plan.beams[0].isocenterPosition_mm, self.overwriteOutsideROI)
        CCCdoseEngineIO.writePlan(self._plan, self._beamDirectory, self.batchSize) 



    def _cleanDir(self, dirPath):
        if os.path.isdir(dirPath):
            shutil.rmtree(dirPath)

    def _startCCC(self, opti=False):
        if len(self._subprocess) > 0:
            raise Exception("CCC already running")

        self._subprocessKilled = False
        logger.info("Start CCC simulation")
        if platform.system() == "Linux":
            for batch in range(self.batchSize):
                if not opti:
                    self._subprocess.append(subprocess.Popen(["sh", 'CCC_simulation_batch{}'.format(batch)], cwd=self._executableDir))
                else:
                    self._subprocess.append(subprocess.Popen(["sh", 'CCC_simulation_opti_batch{}'.format(batch)], cwd=self._executableDir))
            for process in self._subprocess:
                process.wait()

            if self._subprocessKilled:
                self._subprocessKilled = False
                raise Exception('Collapse Cone Convolution subprocess killed by caller.')
            self._subprocess = []


    def _importBeamlets(self):
        beamletDose = CCCdoseEngineIO.readBeamlets(os.path.join(self._ctDirName, 'CT_HeaderFile.txt'), self.outputDir, self.batchSize)
        return beamletDose

    def _importDose(self):
        beamletDose = CCCdoseEngineIO.computeDose(os.path.join(self._ctDirName, 'CT_HeaderFile.txt'), self.outputDir, self.batchSize, self._plan.beamletMUs)
        return beamletDose

    def fromHU2Densities(self, ct : CTImage, overRidingList = None):
        Density = self._ctCalibration._PiecewiseHU2Density__densities
        HU = self._ctCalibration._PiecewiseHU2Density__hu
        ct.imageArray = np.interp(ct.imageArray, HU, Density)
        if overRidingList is not None:
            for overRidingDict in overRidingList:
                ct.imageArray[overRidingDict['Mask'].imageArray.astype(bool) == True] = overRidingDict['Value']
        return ct


    def computeBeamlets(self, ct: CTImage, plan: PhotonPlan, overRidingDict: Optional[Sequence[Union[ROIContour, ROIMask]]] = None) -> SparseBeamlets:
 
        logger.info("Prepare Collapse Cone Convolution Beamlet calculation")
        if self._ct == None:
            self._ct = ct
            self._ct = self.fromHU2Densities(self._ct, overRidingDict) 
            self._plan = plan
        # self._roi = roi

        self._cleanDir(self.outputDir)
        self._cleanDir(self._executableDir)
        self._cleanDir(self._beamDirectory)
        self._writeFilesToSimuDir()
        self.writeExecuteCCCfile()
        self._startCCC()

        beamletDose = self._importBeamlets()

        nbOfBeamlets = beamletDose._sparseBeamlets.shape[1]
        assert(nbOfBeamlets==len(self._plan.beamlets))
        beamletDose.beamletAngles_rad = self._plan.beamletsAngle_rad

        beamletsMU = np.array(plan.beamletMUs)
        if beamletDose.shape[1] != len(beamletsMU):
            print('ERROR: The beamlets imported from the dose engine don\'t have the same number as the beamlets in the plan')
            return
        if plan.planDesign.beamlets is not None:
            beamletDose.beamletWeights = plan.planDesign.beamlets.beamletWeights
        else:
            beamletDose.beamletWeights = beamletsMU
        return beamletDose

    def computeRobustScenarioBeamlets(self, ct: CTImage, plan: PhotonPlan, roi: Optional[Sequence[Union[ROIContour, ROIMask]]] = None, robustMode = "Shift", computeNominal = True) -> SparseBeamlets:
        logger.info("Prepare Collapse Cone Convolution Beamlet calculation")
        self._plan = plan
        self._ct = self.fromHU2Densities(ct, roi)
        self._roi = roi
        self.batchSize = plan.numberOfBeamlets if plan.numberOfBeamlets / self.batchSize < 1 else self.batchSize
        origin = ct.origin
        plan.planDesign.robustness.generateRobustScenarios4Planning()
        scenarios = plan.planDesign.robustness.scenarios
        if computeNominal:
            print('Calculating Nominal Scenario')
            self.ROFolder = 'Nominal'
            plan.planDesign.robustness.nominal.sb = self.computeBeamlets(self._ct, self._plan, roi)   
            plan.planDesign.beamlets = plan.planDesign.robustness.nominal.sb 
        else:
            plan.planDesign.robustness.nominal.sb = plan.planDesign.beamlets 
            
            
        for number, scenario in enumerate(scenarios):
            print('Calculating Scenario {}'.format(number))
            print(scenario)
            self.ROFolder = 'Scenario_{}'.format(number)
            scenario.sb = self.calculateRobustBeamlets(scenario, origin, plan.planDesign.robustness.nominal.sb, mode = robustMode)
        self._ct.origin = origin

    def calculateRobustBeamlets(self, scenario, origin, nominal = None, mode = "Simulation"):
        t0 = time.time()
        self._ct.origin = origin + scenario.sse
        scenario.sre = None if scenario.sre == [0,0,0] else scenario.sre 

        if mode == "Simulation":
            print(origin, self._ct.origin)
            beamletsScenario = self.computeBeamlets(self._ct, self._plan, self._roi)
            beamletsScenario.doseOrigin = origin
        elif mode == "Shift" or scenario.sre != None:
            # kernel = gaussian_kernel_3d(3, scenario.sre)
            scenarioShift_voxel = scenario.sse / self._ct.spacing
            BeamletMatrix = []
            if nominal == None:
                KeyError('To calculate the robust scenarios beamlets in precise mode it is necessary the nominal beamlets')

            nbOfBeamlets = nominal._sparseBeamlets.shape[1]
            assert(nbOfBeamlets==len(self._plan.beamlets))

            BeamletMatrix = shiftBeamlets(nominal._sparseBeamlets, nominal.doseGridSize, scenarioShift_voxel, self._plan.beamletsAngle_rad) ### Implement the convolutions in case of sre in GPU look at shiftBeamlets
            beamletsScenario = SparseBeamlets()
            beamletsScenario.setUnitaryBeamlets(BeamletMatrix)
            beamletsScenario.doseOrigin = nominal.doseOrigin
            beamletsScenario.doseSpacing = nominal.doseSpacing
            beamletsScenario.doseGridSize = nominal.doseGridSize
            beamletsScenario.beamletWeights = nominal.beamletWeights
        else:
            KeyError('The only modes available to calculate the setup scenarios are "Simulation" or "Shift"')
        print('The scenario runned in ',time.time()-t0)
        return beamletsScenario
    
    def computeRobustScenario(self, ct: CTImage, plan: PhotonPlan, roi: Optional[Sequence[Union[ROIContour, ROIMask]]] = None, robustMode = "Shift", computeNominal = True):
        """
        Compute robustness scenario using Collapse Cone Convolution

        Parameters
        ----------


        Returns
        -------
        scenarios:Robustness
            Robustness with nominal and error scenarios
        """
        self._plan = plan
        self._ct = self.fromHU2Densities(ct, roi)
        self._ct = ct
        self._roi = roi
        self.batchSize = plan.numberOfBeamlets if plan.numberOfBeamlets / self.batchSize < 1 else self.batchSize
        origin = ct.origin
        plan.planDesign.robustnessEval.generateRobustScenarios4Planning()
        scenarios = plan.planDesign.robustnessEval.scenarios

        if computeNominal:
            print('Calculating Nominal Scenario')
            self.ROFolder = 'Nominal'
            dose = self.computeDose(self._ct, self._plan, roi)
            #
            plan.planDesign.robustnessEval.nominal.sb = plan.planDesign.robustness.nominal.sb
            plan.planDesign.robustnessEval.nominal.sb.beamletWeights = plan.beamletMUs
            nominal_sb = plan.planDesign.robustnessEval.nominal.sb
            dose = self.calculateDoseArray(nominal_sb, nominal_sb.beamletWeights, 1)
            plan.planDesign.robustnessEval.nominal.dose = dose
        
        for number, scenario in enumerate(scenarios):
            print('Calculating Scenario {}'.format(number))
            print(scenario)
            self.ROFolder = 'Scenario_{}'.format(number)

            scenario.sre = np.array(plan.planDesign.robustnessEval.setupRandomError)
            nominal_scenario = plan.planDesign.robustnessEval.nominal
            plan.planDesign.robustnessEval.scenarios[number].sb = plan.planDesign.robustnessEval.nominal.sb
            plan.planDesign.robustnessEval.scenarios[number].sb.beamletWeights = plan.planDesign.robustnessEval.nominal.sb.beamletWeights
            scenario.dose = self.computeRobustScenarioDose(scenario, origin, nominal_scenario, mode = robustMode)

        self._ct.origin = origin

        return plan.planDesign.robustnessEval

    def calculateDoseArray(self, beamlets, weights, numberOfFractionsPlanned):
        doseArray  = csc_matrix.dot(beamlets._sparseBeamlets, weights) * numberOfFractionsPlanned
        totalDose = np.reshape(doseArray, beamlets._gridSize, order='F')
        totalDose = np.flip(totalDose, 0)
        totalDose = np.flip(totalDose, 1)
        orientation = (1, 0, 0, 0, 1, 0, 0, 0, 1)
        doseImage = DoseImage(imageArray=totalDose, origin=self._ct.origin, spacing=self._ct.spacing,
                        angles=orientation)
        return doseImage

    def _createFolderIfNotExists(self, folder):
        folder = Path(folder)
        if not folder.is_dir():
            os.makedirs(folder)

    def computeDose(self, ct: CTImage, plan: PhotonPlan, overRidingDict: Optional[Sequence[Union[ROIContour, ROIMask]]] = None, Density = False) -> SparseBeamlets:
        logger.info("Prepare Collapse Cone Convolution Beamlet calculation")
        self._ct = ct
        self._plan = plan
        # self._roi = roi

        self._cleanDir(self.outputDir)
        self._cleanDir(self._executableDir)
        self._cleanDir(self._beamDirectory)
        if not Density:
            self._ct = self.fromHU2Densities(self._ct, overRidingDict) 
        self._writeFilesToSimuDir()
        self.writeExecuteCCCfile()
        self._startCCC()


        Dose = self._importDose()
        Dose.imageArray *= self._plan.numberOfFractionsPlanned
        return Dose
  
    def computeRobustScenarioDose(self, scenario, origin, nominal, mode = "Simulation"):
        t0 = time.time()
        self._ct.origin = origin + scenario.sse
        nominal.sb.spacing = self._ct.spacing
        scenario.sre = None if np.all(scenario.sre == [0,0,0]) else scenario.sre

        if mode == "Simulation":
            print(origin, self._ct.origin)
            DoseScenario = self.computeDose(self._ct, self._plan, self._roi)
            DoseScenario.doseOrigin = origin

        elif mode == "Shift" or scenario.sre != None:
            if nominal == None:
                KeyError('To calculate the robust scenarios beamlets in precise mode it is necessary the nominal beamlets')

            nbOfBeamlets = nominal.sb._sparseBeamlets.shape[1]
            assert(nbOfBeamlets==len(self._plan.beamlets))
            DoseScenario = adjustDoseToScenario(scenario, nominal, nominal.sb.spacing, self._plan)

        else:
            KeyError('The only modes available to calculate the setup scenarios are "Simulation" or "Shift"')
        print('The scenario runned in ',time.time()-t0)
        return DoseScenario
