from opentps.core.processing.planOptimization.objectives.dosimetricObjectives._dosimetricObjective import DosimetricObjective
import numpy as np
import scipy.sparse as sp
from scipy import ndimage
import copy

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
    
import logging
logger = logging.getLogger(__name__)

class DFallOff(DosimetricObjective):
    """
    This function defines a distance-dependent maximum dose constraint.

    The maximum permissible dose decreases linearly with
    increasing distance from the given ROI (must be the target).
    It is equal to the prescribed high-dose limit at the target boundary and
    falls to the specified low-dose limit at the designed fall off distance
    from the target.

    The constraint is considered satisfied when all voxel doses within
    the ROI bounded by the target and the fall off contour remain less
    than or equal to their corresponding distance-dependent maximum dose values.

    Attributes
    ----------
    target : image3D
        The target region from which the fall off distance is measured.
    roi : image3D
        The ROI where the DFallOff objective is applied.
    fallOffHighDoseLevel : float
        The maximum dose level at the target boundary.
    fallOffLowDoseLevel : float
        The minimum dose level at the fall off distance from the target.
    fallOffDistance : float
        The distance from the target at which the dose falls to the low-dose level in mm.
    weight : float, optional
        The weight of the objective function, default is 1.
    robust : bool, optional
        If True, the objective function is calculated in a robust way, default is False.
    """
    def __init__(self,target:ROIMask, roi:ROIMask, fallOffHighDoseLevel:float, fallOffLowDoseLevel:float, fallOffDistance:float, weight = 1, robust = False):
        super().__init__(metric='DFallOff',roi=roi,weight=weight,robust=robust)
        if fallOffHighDoseLevel < 0:
            raise ValueError("fallOffHighDoseLevel must be greater or equal to 0 and is currently set to {}".format(fallOffHighDoseLevel))
        self.fallOffHighDoseLevel = fallOffHighDoseLevel
        if fallOffLowDoseLevel < 0:
            raise ValueError("fallOffLowDoseLevel must be greater or equal to 0 and is currently set to {}".format(fallOffLowDoseLevel))
        self.fallOffLowDoseLevel = fallOffLowDoseLevel
        if fallOffDistance < 0:
            raise ValueError("fallOffDistance must be greater or equal to 0 and is currently set to {}".format(fallOffDistance))
        self.fallOffDistance = fallOffDistance

        self.target = target
        self._updateMaskVec(spacing=target.spacing, gridSize=target.gridSize, origin=target.origin)

    def _updateMaskVec(self,spacing: Sequence[float], gridSize: Sequence[int], origin: Sequence[float]):
        targetMask = self.target
        mask = self.roi
        if isinstance(targetMask, ROIContour):
            targetMask = targetMask.getBinaryMask(origin=origin, gridSize=gridSize, spacing=spacing)
        elif isinstance(targetMask, ROIMask):
            targetMask = targetMask
            if not (np.array_equal(targetMask.gridSize, gridSize) and
                    np.allclose(targetMask.origin, origin, atol=0.01) and
                    np.allclose(targetMask.spacing, spacing, atol=0.01)):
                targetMask = resampler3D.resampleImage3D(targetMask, gridSize=gridSize, spacing=spacing, origin=origin)

        euclidDist = ndimage.distance_transform_edt(targetMask.imageArray == 0,sampling=self.target.spacing)  # sampling to express distance in metric units (mm)
        # Check euclid dist size
        masknan = ~np.bool(copy.deepcopy(mask.imageArray))
        masknan[masknan] = 0
        euclidDistROI = euclidDist * masknan

        # Extract voxels within distance in the specified ROI
        voxelsIN = np.logical_and(euclidDistROI > 0, euclidDistROI < self.fallOffDistance)  # ?
        self.maskVec = np.flip(voxelsIN, (0, 1))
        self.maskVec = np.ndarray.flatten(self.maskVec, 'F')
        # get dose rate
        doseRate = (self.fallOffHighDoseLevel - self.fallOffLowDoseLevel) / self.fallOffDistance
        # get reference dose (Dref) as voxel-by-voxel objective
        self.voxelwiseLimitValue = (self.fallOffHighDoseLevel - euclidDistROI * doseRate)  # * voxelsIN  #(self.targetPrescription - euclidDistROI * doseRate) #
        self.voxelwiseLimitValue = np.flip(self.voxelwiseLimitValue, (0, 1))
        # convert 1D vector
        self.voxelwiseLimitValue = np.ndarray.flatten(self.voxelwiseLimitValue, 'F')
        self.voxelwiseLimitValue = self.voxelwiseLimitValue[self.maskVec]

    def _eval(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")

        if self.GPU_acceleration:
            f = cp.mean(cp.maximum(0, dose[self.maskVec_GPU] - self.voxelwiseLimitValue) ** 2)
            self.fValue = f
        else:
            f = np.mean(np.maximum(0, dose[self.maskVec] - self.voxelwiseLimitValue) ** 2)
            self.fValue = f
        return self.fValue

    def _grad(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")

        if self.GPU_acceleration:
            voxelwiseLimitValue_gpu = cp.asarray(self.voxelwiseLimitValue).astype(cp.float32)
            dfdD = 2/cp.sum(self.maskVec_GPU) * cp.maximum(0, dose[self.maskVec_GPU] - voxelwiseLimitValue_gpu)
        else:
            dfdD = 2/np.sum(self.maskVec) * np.maximum(0, dose[self.maskVec] - self.voxelwiseLimitValue)

        if kwargs.get('return_dfdD', False):
            self.gradVector = dfdD
            return self.gradVector
        else:
            dDdx = kwargs.get('dDdx', None)
            if dDdx is None:
                logger.error("Beamlet matrix must be provided")
                raise ValueError("Beamlet matrix must be provided")

        if self.GPU_acceleration:
            dfdx = cpx.scipy.sparse.csc_matrix.dot(dDdx[:, self.maskVec_GPU], dfdD)
            self.gradVector = dfdx
        elif self.MKL_acceleration:
            dfdD = sp.csc_array(dfdD)
            dfdx = sparse_dot_mkl.dot_product_mkl(dDdx[:, self.maskVec], dfdD)
            self.gradVector = dfdx
        else:
            dfdx = sp.csc_matrix.dot(dDdx[:, self.maskVec], dfdD)
            self.gradVector = dfdx
        return self.gradVector