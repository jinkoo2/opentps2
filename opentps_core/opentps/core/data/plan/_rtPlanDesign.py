
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
        self.targetMask: ROIMask = None
        self.calibration: AbstractCTCalibration = None
        self.ct: CTImage = None
        self.beamNames = []
        self.gantryAngles = []
        self.couchAngles = []

        self.objectives = ObjectivesList()
        self.beamlets = []

        self.robustness = Robustness()

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

            
