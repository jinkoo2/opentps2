import os
from opentps.core.data.images import CTImage
from opentps.core.io.scannerReader import readScanner
from opentps.core.processing.doseCalculation.doseCalculationConfig import DoseCalculationConfig
from opentps.core.processing.doseCalculation.photons.cccDoseCalculator import CCCDoseCalculator
from opentps.core.io.sitkIO import exportImageSitk
from opentps.core.data.plan._photonPlan import PhotonPlan
from opentps.core.data.plan._planPhotonBeam import PlanPhotonBeam
from opentps.core.data.plan._planPhotonSegment import PlanPhotonSegment
import numpy as np
from opentps.core.data.images import ROIMask
from opentps.core.data.plan import PhotonPlanDesign
from opentps.core.processing.planEvaluation.robustnessPhotons import Robustness as RobustnessPhotons
from opentps.core.data.plan import ObjectivesList

def create_sphere_in_array(array, radius, center):
    # Get array shape
    shape = array.shape

    # Create indices grid
    indices = np.indices(shape)

    # Calculate distances from the center
    distances = np.sqrt(np.sum((indices - np.array(center).reshape(3, 1, 1, 1))**2, axis=0))

    # Set voxels within the sphere to 1
    array[distances <= radius] = 1

    return array

ctCalibration = readScanner(DoseCalculationConfig().scannerFolder)

size = np.array([300,300,300])
spacing_mm = np.array([2, 2, 2])
phamtom_size = np.array([200,200,200])
Origin = -(size/2 * spacing_mm) 
Origin[1] += 200 

### Create Phantom
Density = 1 ## Water density
img = np.zeros(size)
img[size[0]//2-phamtom_size[0]//2:size[0]//2+phamtom_size[0]//2, size[1]//2-phamtom_size[1]//2:size[1]//2+phamtom_size[1]//2, size[2]//2-phamtom_size[2]//2:size[2]//2+phamtom_size[2]//2] = Density
phantom = CTImage(img,'Phantom',Origin, spacing_mm)
# phantomOutputPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', "phantom.nrrd") 
# exportImageMHD(phantomOutputPath, phantom)

### Create Target Mask
phantomMask = ROIMask('Body')
phantomMask.imageArray = img > 0
mask = 1 
radius = 20
img = np.zeros(size)
array = create_sphere_in_array(img, radius, center = phamtom_size//2)
targetMask = ROIMask('target')
targetMask.imageArray = array 


planInit = PhotonPlanDesign()
planInit.ct = phantom
planInit.targetMask = targetMask
planInit.gantryAngles = np.linspace(0,360,10)
planInit.beamNames = ["Beam{}".format(i) for i in range(len(planInit.gantryAngles))]
planInit.calibration = ctCalibration
planInit.xBeamletSpacing_mm = 5
planInit.yBeamletSpacing_mm = 5
planInit.robustness = RobustnessPhotons()
planInit.robustness.setupSystematicError = [1.6] * 3
planInit.robustness.setupRandomError = 0
planInit.robustness.sseNumberOfSamples = 2
planInit.isocenterPosition_mm = [0,0,0]
plan = planInit.buildPlan() 


ccc = CCCDoseCalculator(batchSize= 15)
ccc.ctCalibration = ctCalibration
plan.planDesign.beamlets = ccc.computeBeamlets(phantom, plan, Density = True)
ccc.computeRobustScenarioBeamlets(phantom, plan, Density = True)

plan.planDesign.objectives = ObjectivesList()
plan.planDesign.objectives.setTarget('target', 70)
plan.planDesign.objectives.addFidObjective(targetMask, 'DMax', 70, 10, robust = True)
plan.planDesign.objectives.addFidObjective(targetMask, 'DMin', 70.5, 10, robust = True)
plan.planDesign.objectives.addFidObjective(phantomMask, 'DMax', 0, 0.1, robust = False)
    
dose_file_nrrd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', "PhantomSSD100cmSAD100cm_1mm_100x100.nrrd") 
exportImageSitk(dose_file_nrrd, 'DMax', 70, 10, robust = True)

