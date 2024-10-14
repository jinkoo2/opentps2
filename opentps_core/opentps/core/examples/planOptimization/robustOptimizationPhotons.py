import os
import logging
import numpy as np
from matplotlib import pyplot as plt
import sys
import scipy as sp
sys.path.append('..')

from opentps.core.data.images import CTImage
from opentps.core.data.images import ROIMask
from opentps.core.data.plan import ObjectivesList
from opentps.core.data.plan._ionPlanDesign import IonPlanDesign
from opentps.core.data import DVH
from opentps.core.data import Patient
from opentps.core.data.plan import FidObjective
from opentps.core.io import mcsquareIO
from opentps.core.io.scannerReader import readScanner
from opentps.core.io.serializedObjectIO import loadRTPlan, saveRTPlan
from opentps.core.processing.doseCalculation.doseCalculationConfig import DoseCalculationConfig
from opentps.core.processing.imageProcessing.resampler3D import resampleImage3DOnImage3D
from opentps.core.processing.planOptimization.planOptimization import IMPTPlanOptimizer
from opentps.core.processing.doseCalculation.photons.cccDoseCalculator import CCCDoseCalculator
from opentps.core.data.plan import PhotonPlanDesign
import copy
from scipy.sparse import csc_matrix
from opentps.core.processing.planEvaluation.robustnessPhotons import Robustness as RobustnessPhotons
from opentps.core.io.dicomIO import writeRTDose


def calculateDoseArray(beamlets, weights, numberOfFractionsPlanned):
    doseArray  = csc_matrix.dot(beamlets._sparseBeamlets, weights) * numberOfFractionsPlanned
    totalDose = np.reshape(doseArray, beamlets._gridSize, order='F')
    totalDose = np.flip(totalDose, 0)
    totalDose = np.flip(totalDose, 1)
    return totalDose

logger = logging.getLogger(__name__)

# Generic example: box of water with squared target
def run(output_path=""):
    if(output_path != ""):
        output_path = output_path
    else:
        output_path = os.path.join(os.getcwd(), 'Photon_Robust_Output_Example')
        if not os.path.exists(output_path):
            os.makedirs(output_path)
    logger.info('Files will be stored in {}'.format(output_path))

    ctCalibration = readScanner(DoseCalculationConfig().scannerFolder)

    patient = Patient()
    patient.name = 'Patient'

    ctSize = 150

    ct = CTImage()
    ct.name = 'CT'
    ct.patient = patient
    
    huAir = -1024.
    huWater = 0
    data = huAir * np.ones((ctSize, ctSize, ctSize))
    data[:, 50:, :] = huWater
    ct.imageArray = data

    roi = ROIMask()
    roi.patient = patient
    roi.name = 'TV'
    roi.color = (255, 0, 0)  # red
    data = np.zeros((ctSize, ctSize, ctSize)).astype(bool)
    data[100:120, 100:120, 100:120] = True
    roi.imageArray = data

    # Design plan
    beamNames = ["Beam1", "Beam2"]
    gantryAngles = [0., 90.]
    couchAngles = [0.,0]

    # Create output folder
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    ## Dose computation from plan
    ccc = CCCDoseCalculator(batchSize= 30)
    ccc.ctCalibration = readScanner(DoseCalculationConfig().scannerFolder)


    # Load / Generate new plan
    plan_file = os.path.join(output_path, "RobustPlan_notCropped.tps")

    if os.path.isfile(plan_file):
        plan = loadRTPlan(plan_file, 'photon')
        logger.info('Plan loaded')
    else:
        planDesign = PhotonPlanDesign()
        planDesign.ct = ct
        planDesign.targetMask = roi
        planDesign.gantryAngles = gantryAngles
        planDesign.beamNames = beamNames
        planDesign.couchAngles = couchAngles
        planDesign.calibration = ctCalibration
        planDesign.xBeamletSpacing_mm = 5
        planDesign.yBeamletSpacing_mm = 5
        # Robustness settings
        planDesign.robustness = RobustnessPhotons()
        planDesign.robustness.setupSystematicError = [1.6] * 3
        # planDesign.robustness.setupRandomError = 1.6
        planDesign.robustness.sseNumberOfSamples = 1

        planDesign.robustness.selectionStrategy = planDesign.robustness.Strategies.REDUCED_SET
        # planDesign.robustness.NumScenarios = 10

        planDesign.targetMargin = max(planDesign.robustness.setupSystematicError)
        planDesign.defineTargetMaskAndPrescription(target = roi, targetPrescription = 20.) # needs to be called prior spot placement
        plan = planDesign.buildPlan()  # Spot placement
        plan.PlanName = "RobustPlan"

        ccc.computeRobustScenarioBeamlets(ct, plan, robustMode='Shift')
        


    saveRTPlan(plan, plan_file, unloadBeamlets=False)
    plan.planDesign.objectives.addFidObjective(roi, FidObjective.Metrics.DMAX, 20.0, 1.0, robust=True)
    plan.planDesign.objectives.addFidObjective(roi, FidObjective.Metrics.DMIN, 20.5, 1.0, robust=True)

    solver = IMPTPlanOptimizer(method='Scipy_L-BFGS-B', plan=plan, maxit=50)
    # Optimize treatment plan
    doseInfluenceMatrix = copy.deepcopy(plan.planDesign.beamlets)
    doseImage, ps = solver.optimize()

    doseImage.imageArray  = calculateDoseArray(doseInfluenceMatrix, plan.beamletMUs, plan.numberOfFractionsPlanned)
    # User input filename
    # writeRTDose(doseImage, output_path, outputFilename="BeamletTotalDose")
    # or default name
    writeRTDose(doseImage, output_path)

    plan_file = os.path.join(output_path, "Plan_Photon_WaterPhantom_cropped_optimized.tps")
    saveRTPlan(plan, plan_file, unloadBeamlets=False)

    # Compute DVH
    target_DVH = DVH(roi, doseImage)
    print('D95 = ' + str(target_DVH.D95) + ' Gy')
    print('D5 = ' + str(target_DVH.D5) + ' Gy')
    print('D5 - D95 =  {} Gy'.format(target_DVH.D5 - target_DVH.D95))

    # center of mass
    roi = resampleImage3DOnImage3D(roi, ct)
    COM_coord = roi.centerOfMass
    COM_index = roi.getVoxelIndexFromPosition(COM_coord)
    Z_coord = COM_index[2]

    img_ct = ct.imageArray[:, :, Z_coord].transpose(1, 0)
    contourTargetMask = roi.getBinaryContourMask()
    img_mask = contourTargetMask.imageArray[:, :, Z_coord].transpose(1, 0)
    img_dose = resampleImage3DOnImage3D(doseImage, ct)
    img_dose = img_dose.imageArray[:, :, Z_coord].transpose(1, 0)

    # Display dose
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].axes.get_xaxis().set_visible(False)
    ax[0].axes.get_yaxis().set_visible(False)
    ax[0].imshow(img_ct, cmap='gray')
    ax[0].imshow(img_mask, alpha=.2, cmap='binary')  # PTV
    dose = ax[0].imshow(img_dose, cmap='jet', alpha=.2)
    plt.colorbar(dose, ax=ax[0])
    ax[1].plot(target_DVH.histogram[0], target_DVH.histogram[1], label=target_DVH.name)
    ax[1].set_xlabel("Dose (Gy)")
    ax[1].set_ylabel("Volume (%)")
    plt.grid(True)
    plt.legend()
    plt.savefig(os.path.join(output_path, 'Dose_RobustOptimizationPhotons.png'))
    # plt.show()

if __name__ == "__main__":
    run()
