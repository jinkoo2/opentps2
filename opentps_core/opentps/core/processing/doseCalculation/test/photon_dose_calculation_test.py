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
from opentps.core.data.plan import PlanDesign, PlanDesignPhotons

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
mask = 1 
radius = 20
img = np.zeros(size)
array = create_sphere_in_array(img, radius, center = phamtom_size//2)
targetMask = ROIMask('target')
targetMask.imageArray = array 

planInit = PlanDesignPhotons()
planInit.ct = ct
planInit.targetMask = targetMask
planInit.gantryAngles = getGantryAngles(config)

plan = PhotonPlan()
beam = PlanPhotonBeam()
plan.SAD_mm = 1000
plan.targetMask 
beam.isocenterPosition_mm = [0,0,0]
segment = PlanPhotonSegment()
segment.xBeamletSpacing_mm = 5
segment.yBeamletSpacing_mm = 5
beamSize = 100 ## Size of a square beam in mm
xRange = np.arange(-beamSize/2,beamSize/2,segment.xBeamletSpacing_mm) + segment.xBeamletSpacing_mm / 2
yRAnge = np.arange(-beamSize/2,beamSize/2,segment.xBeamletSpacing_mm) + segment.yBeamletSpacing_mm / 2
# segment.appendBeamlet(0, 0, 1)

for x in xRange:
    for y in yRAnge:
        segment.appendBeamlet(x, y, 1)
beam.appendBeamSegment(segment)
plan.appendBeam(beam)

ccc = CCCDoseCalculator(batchSize= 20)
ccc.ctCalibration = ctCalibration
dose = ccc.calculateDose(phantom, plan, Density = True)

dose_file_nrrd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', "PhantomSSD100cmSAD100cm_1mm_100x100.nrrd") 
exportImageSitk(dose_file_nrrd, dose)

