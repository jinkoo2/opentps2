from opentps.core.processing.planOptimization.objectives.baseFunction import BaseFunc
import logging
import numpy as np
import scipy.sparse as sp
from opentps.core.data.images._image3D import Image3D
from opentps.core.data.images._roiMask import ROIMask
from opentps.core.processing.imageProcessing import resampler3D
from opentps.core.data._roiContour import ROIContour
from typing import Optional, Sequence, Union, Iterable
try:
    import sparse_dot_mkl
    sdm_available = True
except:
    sdm_available = False

try:
    import mkl as mkl
    mkl_available = True
except:
    mkl_available = False

try:
    import cupy as cp
    import cupyx as cpx
    cupy_available = True
except:
    cupy_available = False

logger = logging.getLogger(__name__)
from enum import Enum

class DosimetricObjective(BaseFunc):
    """
    Base class for dosimetric objectives.
    This class serves as a foundation for specific dosimetric objectives used in treatment plan optimization.
    It handles common functionalities such as dose retrieval, ROI mask management, and GPU acceleration support.
    It also handles the update of the dose at each iteration because of the change of the beamlet weights.

    Attributes:
        metric : str
            The dosimetric metric to be optimized (e.g., 'DMAX', 'DMIN',...).
        roi : Union[ROIContour, ROIMask]
            The region of interest (ROI) for which the dosimetric objective is defined.
        dose : Optional[Union[np.ndarray, cp.ndarray]]  (default: None)
            The dose distribution associated with the current beamlet weights.
        maskVec : Optional[np.ndarray] (default: None)
            A flattened binary mask vector representing the ROI.
        maskVec_GPU : Optional[cp.ndarray] (default: None)
            A GPU-accelerated version of the mask vector, if GPU acceleration is enabled.
        kind : str (default: "Soft")
            The type of dosimetric objective. Will be deprecated in future versions.
    """
    def __init__(self,metric,roi,robust = False,weight = 1):
        super().__init__(robust = robust, weight = weight)
        self.metric = metric
        self.roi = roi
        self.dose = None
        self.maskVec = None
        self.maskVec_GPU = None
        self.kind = "Soft"

    class Metrics(Enum):
        DMIN = 'DMin'
        DMAX = 'DMax'
        DMAXMEAN = 'DMaxMean'
        DMINMEAN = 'DMinMean'
        DUNIFORM = 'DUniform'
        DVHMIN = 'DVHMin'
        DVHMAX = 'DVHMax'
        DFALLOFF = 'DFallOff'
        EUDMIN = 'EUDMin'
        EUDMAX = 'EUDMax'
        EUDUNIFORM = 'EUDUniform'

    def _eval(self, x, **kwargs):
        logger.error("This dosimetric objective does not have an evaluation function")
        return None

    def _grad(self, x, **kwargs):
        logger.error("This dosimetric objective does not have a gradient function")
        return None

    def _calcInverseDVH(self,volume,dose):
        if self.GPU_acceleration:
            sorted_dose = cp.sort(dose.flatten())
            index = cp.rint((1 - volume) * len(sorted_dose)).astype(cp.int32)
        else:
            sorted_dose = np.sort(dose.flatten())
            index = int((1 - volume) * len(sorted_dose))
        return sorted_dose[index]

    def _updateMaskVec(self, spacing: Sequence[float], gridSize: Sequence[int], origin: Sequence[float]):

        if isinstance(self.roi, ROIContour):
            mask = self.roi.getBinaryMask(origin=origin, gridSize=gridSize, spacing=spacing)
        elif isinstance(self.roi, ROIMask):
            mask = self.roi
            if not (np.array_equal(mask.gridSize, gridSize) and
                    np.allclose(mask.origin, origin, atol=0.01) and
                    np.allclose(mask.spacing, spacing, atol=0.01)):
                mask = resampler3D.resampleImage3D(self.roi, gridSize=gridSize, spacing=spacing, origin=origin)
        else:
            raise Exception(self.roi.__class__.__name__ + ' is not a supported class for roi')

        self.maskVec = np.flip(mask.imageArray, (0, 1))
        self.maskVec = np.ndarray.flatten(self.maskVec, 'F').astype('bool')

    def _loadMaskVecToGPU(self):
        self.maskVec_GPU = cp.asarray(self.maskVec)

    @property
    def roiName(self) -> str:
        return self.roi.name


