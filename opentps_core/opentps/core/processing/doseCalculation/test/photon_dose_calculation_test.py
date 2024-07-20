import os
from opentps.core.data.images._ctImage import CTImage
from opentps.core.io.scannerReader import readScanner
from opentps.core.processing.doseCalculation.doseCalculationConfig import DoseCalculationConfig
from opentps.core.processing.doseCalculation.photons.cccDoseCalculator import CCCDoseCalculator
from opentps.core.io.sitkIO import exportImageSitk
import numpy as np
from opentps.core.data.images import ROIMask
from opentps.core.data.plan import PhotonPlanDesign
from opentps.core.processing.planEvaluation.robustnessPhotons import Robustness as RobustnessPhotons
from opentps.core.data.plan import ObjectivesList
from opentps.core.processing.planOptimization.planOptimization import IMPTPlanOptimizer
import pickle
import glob
import copy


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

def joinDir(path1,path2=""):
    dir = os.path.join(path1,path2)
    if not os.path.isdir(dir):
        os.makedirs(dir)
    return dir

def saveDIM(plan, output_path, scenarioMode, elementToSave):
    MatrixPath = joinDir(output_path,'SparseMatrices_'+scenarioMode)
    if str(elementToSave) in ['Nominal', 'all', 'True']:
        if plan.planDesign.robustness is None:
            plan.planDesign.beamlets.storeOnFS(os.path.join(MatrixPath,'SM_nominal.pkl'))
            scenarios = []
        else:
            plan.planDesign.robustness.nominal.sb.storeOnFS(os.path.join(MatrixPath,'SM_nominal.pkl'))
            scenarios = plan.planDesign.robustness.scenarios
    if str(elementToSave) in ['Scenarios', 'all', 'True']:    
        for i, scenario in enumerate(scenarios):
            dose_file = os.path.join(MatrixPath,'SM_scenario_{}_{}.pkl'.format(np.round(scenario.sse,1), scenario.sre))
            scenario.sb.storeOnFS(dose_file) 

def loadBeamlets(file_path):
    from opentps.core.data._sparseBeamlets import SparseBeamlets
    return loadData(file_path, SparseBeamlets)

def loadData(file_path, cls):
    with open(file_path, 'rb') as fid:
        tmp = pickle.load(fid)
    data = cls()
    data.__dict__.update(tmp)
    return data

def loadDIM(plan, path):
    if os.path.isdir(path):
        if os.path.isfile(os.path.join(path, 'SM_nominal.pkl')):
            plan.planDesign.beamlets = loadBeamlets(os.path.join(os.path.join(path, 'SM_nominal.pkl')))
            plan.planDesign.robustness.setNominal(plan.planDesign.beamlets)
        scenarioPaths = glob.glob(os.path.join(path,'*scenario*')) 
        if len(scenarioPaths)>0:
            for scenario in scenarioPaths:
                fileName = os.path.basename(scenario)
                if fileName.split('_')[2]!='None' and '[' in fileName and ']' in fileName:
                    split = fileName.split('[')[1].split(']')[0].split(' ')
                    sse = [float(s) for s in split if len(s) > 0 and (s[0].isnumeric() or s[0]=='-')]
                    assert(len(sse) == 3)
                else:
                    sse = None
                if fileName.split('_')[3].split('.')[0]!='None':
                    sre = float(fileName.split('_')[3])
                else:
                    sre = None       
                print(f'loading scenario... sse = {sse} sre = {sre}')         
                plan.planDesign.robustness.addScenario(loadBeamlets(scenario), sse, sre)

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
exportImageSitk(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'phantom.nrrd'), phantom)

### Create Target Mask
phantomMask = ROIMask('Body')
phantomMask.imageArray = img > 0
mask = 1 
radius = 20
img = np.zeros(size)
array = create_sphere_in_array(img, radius, center = size//2)
targetMask = ROIMask('target')
targetMask.imageArray = array 
targetMask.spacing = spacing_mm
targetMask.origin = Origin
exportImageSitk(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'mask.nrrd'), targetMask, mask=True)


planInit = PhotonPlanDesign()
planInit.ct = phantom
planInit.targetMask = targetMask
planInit.gantryAngles = np.linspace(0,360,10)
planInit.couchAngles = np.zeros(len(planInit.gantryAngles))
planInit.beamNames = ["Beam{}".format(i) for i in range(len(planInit.gantryAngles))]
planInit.calibration = ctCalibration
planInit.xBeamletSpacing_mm = 5
planInit.yBeamletSpacing_mm = 5
planInit.robustness = RobustnessPhotons()
planInit.robustness.setupSystematicError = [1.6] * 3
planInit.robustness.setupRandomError = 0
planInit.robustness.sseNumberOfSamples = 1
plan = planInit.buildPlan() 

DIM_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'SparseMatrices_MinMax')
if not os.path.isdir(DIM_path) or len(os.listdir(DIM_path)) == 0:
    ccc = CCCDoseCalculator(batchSize= 30)
    ccc.ctCalibration = ctCalibration
    plan.planDesign.beamlets = ccc.computeBeamlets(phantom, plan)
    # plan.planDesign.beamlets = loadBeamlets(os.path.join(os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'SparseMatrices_MinMax','SM_nominal.pkl'))))
    plan.planDesign.robustness.nominal.sb = plan.planDesign.beamlets
    # plan.planDesign.beamlets = ccc.computeBeamlets(phantom, plan)
    ccc.computeRobustScenarioBeamlets(phantom, plan, computeNominal = False)
    saveDIM(plan, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results'), 'MinMax', 'all')
else:
    loadDIM(plan, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'SparseMatrices_MinMax'))

doseInfluenceMatrix = copy.deepcopy(plan.planDesign.beamlets)

plan.planDesign.objectives = ObjectivesList()
plan.planDesign.objectives.setTarget('target', 70)
plan.planDesign.objectives.addFidObjective(targetMask, 'DMax', 70, 10, robust = True)
plan.planDesign.objectives.addFidObjective(targetMask, 'DMin', 70.5, 10, robust = True)
plan.planDesign.objectives.addFidObjective(phantomMask, 'DMax', 0, 0.1, robust = False)

solver = IMPTPlanOptimizer(method = 'Scipy-LBFGS', plan=plan,
                            maxit = 200,
                            ftol = 1e-6,
                            gtol = 1e-6)
totalDose, cost = solver.optimize()

# doseInfluenceMatrix.beamletWeights = plan.planDesign.beamlets.beamletWeights
# dose = doseInfluenceMatrix.toDoseImage()
dose_file_nrrd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', "dose.nrrd") 
exportImageSitk(dose_file_nrrd, totalDose)

