import os
from opentps.core.data.images._ctImage import CTImage
from opentps.core.io.scannerReader import readScanner
from opentps.core.processing.doseCalculation.doseCalculationConfig import DoseCalculationConfig
from opentps.core.processing.doseCalculation.photons.cccDoseCalculator import CCCDoseCalculator
from opentps.core.io.sitkIO import exportImageSitk
import numpy as np
from opentps.core.data.images import ROIMask
import logging
from opentps.core.data.plan._photonPlan import PhotonPlan
from opentps.core.data.plan._planPhotonBeam import PlanPhotonBeam
from opentps.core.data.plan._planPhotonSegment import PlanPhotonSegment
from opentps.core.io.serializedObjectIO import loadRTPlan, saveRTPlan
from pathlib import Path
from opentps.core.processing.imageProcessing.resampler3D import resampleImage3DOnImage3D
from opentps.core.data._dvh import DVH
import matplotlib.pyplot as plt
from opentps.core.data import Patient

import os
import datetime
import logging
import copy
import numpy as np

from opentps.core.io.dicomIO import writeRTDose, readDicomDose
from matplotlib import pyplot as plt
from opentps.core.data.images import CTImage
from opentps.core.data.images import ROIMask
from opentps.core.data.plan._ionPlanDesign import IonPlanDesign
from opentps.core.data import Patient
from opentps.core.io import mcsquareIO
from opentps.core.io.scannerReader import readScanner
from opentps.core.io.serializedObjectIO import saveRTPlan, loadRTPlan
from opentps.core.processing.doseCalculation.doseCalculationConfig import DoseCalculationConfig
from opentps.core.processing.planEvaluation.robustnessPhotons import RobustEvaluation
from opentps.core.processing.planEvaluation.robustnessPhotons import Robustness as RobustnessPhotons
from opentps.core.processing.doseCalculation.photons.cccDoseCalculator import CCCDoseCalculator
from opentps.core.data.plan import PhotonPlanDesign
from scipy.sparse import csc_matrix

ctCalibration = readScanner(DoseCalculationConfig().scannerFolder)
logger = logging.getLogger(__name__)


def run(output_path=""):
    if(output_path != ""):
        output_path = output_path
    else:
        output_path = os.path.join(os.getcwd(), 'Output_Example','PhotonDoseCalculation')
        if not os.path.exists(output_path):
            os.makedirs(output_path)

    logger.info('Files will be stored in {}'.format(output_path))

   
    # Load plan
    plan = loadRTPlan('/home/colin/opentps/Photon_Robust_Output_Example/Plan_Photon_WaterPhantom_cropped_optimized.tps')
    print(plan)
    plan.beamNames = ["Beam1", "Beam2"]
    plan.gantryAngles = [0., 90.]
    plan.couchAngles = [0.,0]
    
    
    ## Dose computation from plan
    ccc = CCCDoseCalculator(batchSize= 30)
    ccc.ctCalibration = readScanner(DoseCalculationConfig().scannerFolder)

    patient = Patient()
    patient.name = 'Patient'

    ctSize = 150

    ct = CTImage()
    ct.name = 'CT'
    ct.patient = patient


    huAir = -1024.
    huWater = ctCalibration.convertRSP2HU(1.)
    data = huAir * np.ones((ctSize, ctSize, ctSize))
    data[:, 50:, :] = huWater
    ct.imageArray = data

    roi = ROIMask()
    roi.patient = patient
    roi.name = 'TV'
    roi.color = (255, 0, 0) # red
    data = np.zeros((ctSize, ctSize, ctSize)).astype(bool)
    data[100:120, 100:120, 100:120] = True
    roi.imageArray = data
    
    doseImage = ccc.computeDose(ct, plan)
    
    # DVH
    dvh = DVH(target, doseImage)
    print("D95",dvh._D95)
    print("D5",dvh._D5)
    print("Dmax",dvh._Dmax)
    print("Dmin",dvh._Dmin)
    
    # Plot dose
    target = resampleImage3DOnImage3D(target, ct)
    COM_coord = target.centerOfMass
    COM_index = target.getVoxelIndexFromPosition(COM_coord)
    Z_coord = COM_index[2]

    img_ct = ct.imageArray[:, :, Z_coord].transpose(1, 0)
    contourTargetMask = target.getBinaryContourMask()
    img_mask = contourTargetMask.imageArray[:, :, Z_coord].transpose(1, 0)
    img_dose = resampleImage3DOnImage3D(doseImage, ct)
    img_dose = img_dose.imageArray[:, :, Z_coord].transpose(1, 0)

    # Display dose
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    ax[0].axes.get_xaxis().set_visible(False)
    ax[0].axes.get_yaxis().set_visible(False)
    ax[0].imshow(img_ct, cmap='gray')
    ax[0].imshow(img_mask, alpha=.2, cmap='binary')  # PTV
    dose = ax[0].imshow(img_dose, cmap='jet', alpha=.2)
    plt.colorbar(dose, ax=ax[0])
    ax[1].plot(dvh.histogram[0], dvh.histogram[1], label=dvh.name)
    ax[1].set_xlabel("Dose (Gy)")
    ax[1].set_ylabel("Volume (%)")
    ax[1].grid(True)
    ax[1].legend()
    plt.savefig(os.path.join(output_path, 'dose.png'))    
    # plt.show()    
    
    
if __name__ == "__main__":
    run()