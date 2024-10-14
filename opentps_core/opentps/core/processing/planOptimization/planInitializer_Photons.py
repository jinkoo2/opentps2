import logging
import math

import numpy as np
from opentps.core.data.CTCalibrations._abstractCTCalibration import AbstractCTCalibration
from opentps.core.data.images._ctImage import CTImage
from opentps.core.data.images._roiMask import ROIMask
from opentps.core.data.plan._planIonLayer import PlanIonLayer
from opentps.core.data.plan._planPhotonBeam import PlanPhotonBeam
from opentps.core.processing.C_libraries.libRayTracing_wrapper import transport_spots_to_target, \
    transport_spots_inside_target
from opentps.core.processing.rangeEnergy import energyToRange, rangeToEnergy
from opentps.core.data.plan._photonPlan import PhotonPlan
from opentps.core.data.plan._planPhotonBeam import PlanPhotonBeam

logger = logging.getLogger(__name__)

class BeamInitializer:
    """
    This class is used to initialize a photon beam.

    Attributes
    ----------
    calibration : AbstractCTCalibration
        The CT calibration used to convert the CT image to RSP image.
    targetMargin : float
        The margin around the target in mm.
    beam : PlanPhotonBeam
        The beam to initialize.
    """
    def __init__(self):
        self.targetMargin = 0.
        self.beam: PlanPhotonBeam = None

        self.calibration: AbstractCTCalibration = None

    def initializeBeam(self):
        """
        Initialize the beam with beamlets.
        """
        # generate hexagonal spot grid around isocenter
        beamletGrid = self._definePencilBeamGridAroundIsocenter()
        numBeamlets = len(beamletGrid["x"])

        # compute direction vector
        u, v, w = 1e-10, 1.0, 1e-10  # BEV to 3D coordinates
        [u, v, w] = self._rotateVector([u, v, w], math.radians(self.beam.gantryAngle_degree), 'z')  # rotation for gantry angle
        [u, v, w] = self._rotateVector([u, v, w], math.radians(self.beam.couchAngle_degree), 'y')  # rotation for couch angle

        # prepare raytracing: translate initial positions at the CT image border
        for s in range(numBeamlets):
            translation = np.array([1.0, 1.0, 1.0])
            translation[0] = (beamletGrid["x"][s] - self.imgBordersX[int(u < 0)]) / u
            translation[1] = (beamletGrid["y"][s] - self.imgBordersY[int(v < 0)]) / v
            translation[2] = (beamletGrid["z"][s] - self.imgBordersZ[int(w < 0)]) / w
            translation = translation.min()
            beamletGrid["x"][s] = beamletGrid["x"][s] - translation * u
            beamletGrid["y"][s] = beamletGrid["y"][s] - translation * v
            beamletGrid["z"][s] = beamletGrid["z"][s] - translation * w

        # This function works for protons that is why it requires the RSP. However, if we just take the x,y dimension it is ok for photons. 
        transport_spots_to_target(self.rspImage, self.targetMask, beamletGrid, [u, v, w])

        # remove spots that didn't reach the target
        minWET = 9999999
        for s in range(numBeamlets - 1, -1, -1):
            if beamletGrid["WET"][s] < 0:
                beamletGrid["BEVx"].pop(s)
                beamletGrid["BEVy"].pop(s)
                beamletGrid["x"].pop(s)
                beamletGrid["y"].pop(s)
                beamletGrid["z"].pop(s)
                beamletGrid["WET"].pop(s)
            else:
                if beamletGrid["WET"][s] < minWET: minWET = beamletGrid["WET"][s]

        # process valid spots
        numBeamlets = len(beamletGrid["x"])
        for n in range(numBeamlets):
            self.beam.appendBeamlet(beamletGrid["BEVx"][n], beamletGrid["BEVy"][n], 1)

    def _definePencilBeamGridAroundIsocenter(self):
        FOV = 400  # max field size on IBA P+ is 30x40 cm
        numSpotX = math.ceil(FOV / self.beam.xBeamletSpacing_mm)
        numSpotY = math.ceil(FOV / self.beam.yBeamletSpacing_mm)

        beamletGrid = {"BEVx": [], "BEVy": [], "x": [], "y": [], "z": [], "WET": []}

        for i in range(numSpotX):
            for j in range(numSpotY):
                # coordinates in Beam-eye-view
                beamletGrid["BEVx"].append((i - round(numSpotX / 2)) * self.beam.xBeamletSpacing_mm)
                beamletGrid["BEVy"].append((j - round(numSpotY / 2)) * self.beam.yBeamletSpacing_mm)

                # 3D coordinates
                x, y, z = beamletGrid["BEVx"][-1], 0, beamletGrid["BEVy"][-1]

                # rotation for gantry angle (around Z axis)
                [x, y, z] = self._rotateVector([x, y, z], math.radians(self.beam.gantryAngle_degree), 'z')

                # rotation for couch angle (around Y axis)
                [x, y, z] = self._rotateVector([x, y, z], math.radians(self.beam.couchAngle_degree), 'y')

                # Dicom CT coordinates
                beamletGrid["x"].append(x + self.beam.isocenterPosition_mm[0])
                beamletGrid["y"].append(y + self.beam.isocenterPosition_mm[1])
                beamletGrid["z"].append(z + self.beam.isocenterPosition_mm[2])
        return beamletGrid

    def _rotateVector(self, vec, angle, axis):
        if axis == 'x':
            x = vec[0]
            y = vec[1] * math.cos(angle) - vec[2] * math.sin(angle)
            z = vec[1] * math.sin(angle) + vec[2] * math.cos(angle)
        elif axis == 'y':
            x = vec[0] * math.cos(angle) + vec[2] * math.sin(angle)
            y = vec[1]
            z = -vec[0] * math.sin(angle) + vec[2] * math.cos(angle)
        elif axis == 'z':
            x = vec[0] * math.cos(angle) - vec[1] * math.sin(angle)
            y = vec[0] * math.sin(angle) + vec[1] * math.cos(angle)
            z = vec[2]

        return [x, y, z]


class PhotonPlanInitializer:
    def __init__(self):
        self.ctCalibration: AbstractCTCalibration = None
        self.ct: CTImage = None
        self.plan: PhotonPlan = None
        self.targetMask: ROIMask = None

        self._beamInitializer = BeamInitializer()

    def placeBeamlets(self, targetMargin: float = 0.):
        self._beamInitializer.calibration = self.ctCalibration
        self._beamInitializer.targetMargin = targetMargin

        from opentps.core.data.images._rspImage import RSPImage
        logger.info('Target is dilated using a margin of {} mm. This process might take some time.'.format(targetMargin))
        roiDilated = ROIMask.fromImage3D(self.targetMask, patient=None)
        roiDilated.dilateMask(radius=targetMargin)
        logger.info('Dilation done.')
        self._beamInitializer.targetMask = roiDilated

        rspImage = RSPImage.fromCT(self.ct, self.ctCalibration, energy=100.) ######################################## Review!!!!!!!!!!!!!!!!!!!
        rspImage.patient = None
        self._beamInitializer.rspImage = rspImage

        imgBordersX = [rspImage.origin[0], rspImage.origin[0] + rspImage.gridSize[0] * rspImage.spacing[0]]
        imgBordersY = [rspImage.origin[1], rspImage.origin[1] + rspImage.gridSize[1] * rspImage.spacing[1]]
        imgBordersZ = [rspImage.origin[2], rspImage.origin[2] + rspImage.gridSize[2] * rspImage.spacing[2]]

        self._beamInitializer.imgBordersX = imgBordersX
        self._beamInitializer.imgBordersY = imgBordersY
        self._beamInitializer.imgBordersZ = imgBordersZ

        for beam in self.plan:
            beam.removeBeamSegment(beam.beamSegments)
            if beam.beamType == 'Static':
                beam.createBeamSegment()                
                self._beamInitializer.beam = beam[0]
                self._beamInitializer.initializeBeam()
