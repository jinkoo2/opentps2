
__all__ = ['RTPlanDesign']

import logging
import time
from typing import Optional, Sequence, Union

import numpy as np
import pydicom

from opentps.core.data.CTCalibrations._abstractCTCalibration import AbstractCTCalibration
from opentps.core.data.images._ctImage import CTImage
from opentps.core.data.images._roiMask import ROIMask
from opentps.core.data.plan._rangeShifter import RangeShifter
from opentps.core.processing.imageProcessing import resampler3D
from opentps.core.data._patientData import PatientData
from opentps.core.data.plan._objectivesList import ObjectivesList
from opentps.core.processing.planEvaluation.robustnessEvaluation import Robustness
from opentps.core.processing.planOptimization.planInitializer import PlanInitializer

logger = logging.getLogger(__name__)


class RTPlanDesign(PatientData):
    """
    This class is used to store the plan design. It inherits from PatientData.

    Attributes
    ----------
    targetMargin: float (default: 5.0)
        margin around the target in mm
    targetMask: ROIMask
        mask of the target
    calibration: AbstractCTCalibration
        calibration of the CT for stopping power conversion
    ct: CTImage (default: None)
        CT image
    beamNames: list of str
        list of beam names
    gantryAngles: list of float
        list of gantry angles
    couchAngles: list of float
        list of couch angles
    objectives: ObjectivesList
        list of objectives
    beamlets: list of Beamlet
        list of beamlets
    robustness: Robustness
        robustness evaluation
    """
    def __init__(self):
        super().__init__()

        self.targetMargin = 5.0
        self._scoringVoxelSpacing = None
        self._scoringGridSize = None
        self._scoringOrigin = None
        self.targetMask: ROIMask = None
        self.calibration: AbstractCTCalibration = None
        self.ct: CTImage = None
        self.beamNames = []
        self.gantryAngles = []
        self.couchAngles = []
        self._scoringVoxelSpacing = None

        self.objectives = ObjectivesList()
        self.beamlets = None

        self.robustness = Robustness()


    @property
    def scoringVoxelSpacing(self) -> Sequence[float]:
        if self._scoringVoxelSpacing is not None:
            return self._scoringVoxelSpacing
        else:
            return self.ct.spacing

    @scoringVoxelSpacing.setter
    def scoringVoxelSpacing(self, spacing: Union[float, Sequence[float]]):
        if np.isscalar(spacing):
            self._scoringVoxelSpacing = np.array([spacing, spacing, spacing])
        else:
            self._scoringVoxelSpacing = np.array(spacing)

    @property
    def scoringGridSize(self):
        if self._scoringGridSize is not None:
            return self._scoringGridSize
        else:
            return self.ct.gridSize
    
    @scoringGridSize.setter
    def scoringGridSize(self, gridSize: Sequence[float]):
        self._gridSize = gridSize

    @property
    def scoringOrigin(self):
        if self._scoringOrigin is not None:
            return self._scoringOrigin
        else:
            return self.ct.origin
        
    @scoringOrigin.setter
    def scoringOrigin(self, origin):
        self._scoringOrigin = origin
        
    def defineTargetMaskAndPrescription(self):
        """
        Defines the target mask and the prescription
        """
        from opentps.core.data._roiContour import ROIContour

        targetMask = None
        for objective in self.objectives.fidObjList:
            if objective.metric == objective.Metrics.DMIN:
                roi = objective.roi

                self.objectives.setTarget(objective.roiName, objective.limitValue)

                if isinstance(roi, ROIContour):
                    mask = roi.getBinaryMask(origin=self.ct.origin, gridSize=self.ct.gridSize,
                                             spacing=self.ct.spacing)
                elif isinstance(roi, ROIMask):
                    mask = resampler3D.resampleImage3D(roi, origin=self.ct.origin,
                                                       gridSize=self.ct.gridSize,
                                                       spacing=self.ct.spacing)
                else:
                    raise Exception(roi.__class__.__name__ + ' is not a supported class for roi')

                if targetMask is None:
                    targetMask = mask
                else:
                    targetMask.imageArray = np.logical_or(targetMask.imageArray, mask.imageArray)

        if targetMask is None:
            raise Exception('Could not find a target volume in dose fidelity objectives')

        self.targetMask = targetMask


    def buildPlan(self):
        """
        Builds a plan from the plan design
        """
        pass

    def createBeams(self):
        """
        Creates the beams of the plan

        """
        pass

    def initializeBeams(self):
        """
        Initializes the beams of the plan
        """
        pass
    
    def setScoringParameters(self, scoringGridSize:Optional[Sequence[int]]=None, scoringSpacing:Optional[Sequence[float]]=None,
                                scoringOrigin:Optional[Sequence[int]]=None, adapt_gridSize_to_new_spacing=False):
        """
        Sets the scoring parameters

        Parameters
        ----------
        scoringGridSize: Sequence[int]
            scoring grid size
        scoringSpacing: Sequence[float]
            scoring spacing
        scoringOrigin: Sequence[float]
            scoring origin
        adapt_gridSize_to_new_spacing: bool
            If True, automatically adapt the gridSize to the new spacing
        """
        if adapt_gridSize_to_new_spacing and scoringGridSize is not None:
            raise ValueError('Cannot adapt gridSize to new spacing if scoringGridSize provided.')
        
        if scoringSpacing is not None: self.scoringVoxelSpacing = scoringSpacing
        if scoringGridSize is not None: self.scoringGridSize = scoringGridSize
        if scoringOrigin is not None: self.scoringOrigin = scoringOrigin
        
        if adapt_gridSize_to_new_spacing:
            self.scoringGridSize = np.floor(self.ct.gridSize*self.ct.spacing/self.scoringVoxelSpacing).astype(int)

        for objective in self.objectives.fidObjList:
            objective._updateMaskVec(spacing=self.scoringVoxelSpacing, gridSize=self.scoringGridSize, origin=self.scoringOrigin)
            
