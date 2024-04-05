
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
from opentps.core.data.plan._rtPlanDesign import RTPlanDesign

logger = logging.getLogger(__name__)


class IonPlanDesign(RTPlanDesign):
    """
    This class is used to store the plan design. It inherits from PatientData.

    Attributes
    ----------
    spotSpacing: float (default: 5.0)
        spacing between spots in mm
    layerSpacing: float (default: 5.0)
        spacing between layers in mm
    scoringVoxelSpacing: float or list of float
        spacing of the scoring grid in mm
    proximalLayers: int (default: 1)
        number of proximal layers
    distalLayers: int (default: 1)
        number of distal layers
    layersToSpacingAlignment: bool (default: False)
        if True, the spacing between layers is aligned with the scoring grid
    rangeShifters: list of RangeShifter
        list of range shifters
    beamletsLET: list of Beamlet
        list of beamlets with LET
    """
    def __init__(self):
        super().__init__()

        self.spotSpacing = 5.0
        self.layerSpacing = 5.0
        self._scoringVoxelSpacing = None
        self.proximalLayers = 1
        self.distalLayers = 1
        self.layersToSpacingAlignment = False
        self.rangeShifters: Sequence[RangeShifter] = []

        self.beamletsLET = []

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
        if self._scoringVoxelSpacing is not None:
            return np.floor(self.ct.gridSize*self.ct.spacing/self.scoringVoxelSpacing).astype(int)
        else:
            return self.ct.gridSize

    def buildPlan(self):
        """
        Builds a plan from the plan design

        Returns
        --------
        RTPlan
            plan
        """
        start = time.time()
        # Spot placement
        from opentps.core.data.plan import IonPlan
        plan = IonPlan("NewPlan")
        plan.SOPInstanceUID = pydicom.uid.generate_uid()
        plan.seriesInstanceUID = plan.SOPInstanceUID + ".1"
        plan.modality = "Ion therapy"
        plan.radiationType = "Proton"
        plan.scanMode = "MODULATED"
        plan.treatmentMachineName = "Unknown"
        logger.info('Building plan ...')
        self.createBeams(plan)
        self.initializeBeams(plan)
        plan.planDesign = self
        for beam in plan.beams:
            beam.reorderLayers('decreasing')

        logger.info("New plan created in {} sec".format(time.time() - start))
        logger.info("Number of spots: {}".format(plan.numberOfSpots))

        return plan

    def createBeams(self, plan):
        """
        Creates the beams of the plan

        Parameters
        ----------
        plan: RTPlan
            plan
        """
        for beam in plan:
            plan.removeBeam(beam)

        from opentps.core.data.plan import PlanIonBeam
        for i, gantryAngle in enumerate(self.gantryAngles):
            beam = PlanIonBeam()
            beam.gantryAngle = gantryAngle
            beam.couchAngle = self.couchAngles[i]
            beam.isocenterPosition = self.targetMask.centerOfMass
            beam.id = i
            if self.beamNames:
                beam.name = self.beamNames[i]
            else:
                beam.name = 'B' + str(i)
            if self.rangeShifters and self.rangeShifters[i]:
                beam.rangeShifter = self.rangeShifters[i]

            plan.appendBeam(beam)

    def initializeBeams(self, plan):
        """
        Initializes the beams of the plan

        Parameters
        ----------
        plan: RTPlan
            plan
        """
        initializer = PlanInitializer()
        initializer.ctCalibration = self.calibration
        initializer.ct = self.ct
        initializer.plan = plan
        initializer.targetMask = self.targetMask
        initializer.placeSpots(self.spotSpacing, self.layerSpacing, self.targetMargin, self.layersToSpacingAlignment,
                               self.proximalLayers, self.distalLayers)


    def setScoringParameters(self, scoringGridSize:Optional[Sequence[int]]=None, scoringSpacing:Optional[Sequence[float]]=None):
        """
        Sets the scoring parameters

        Parameters
        ----------
        scoringGridSize: Sequence[int]
            scoring grid size
        scoringSpacing: Sequence[float]
            scoring spacing
        """
        if scoringSpacing is None and scoringGridSize is not None:
            self.scoringVoxelSpacing = self.ct.spacing*self.ct.gridSize/scoringGridSize
        if scoringSpacing is not None and scoringGridSize is None:
            self.scoringVoxelSpacing = scoringSpacing
        if scoringSpacing is not None and scoringGridSize is not None:
            raise Exception('Cannot set both scoring spacing and grid size at the same time.')
        # scoringSpacing and scoringGridSize are None --> defaults to CT spacing and size


        for objective in self.objectives.fidObjList:
            objective._updateMaskVec(spacing=self.scoringVoxelSpacing, gridSize=self.scoringGridSize, origin=self.ct.origin)
            
