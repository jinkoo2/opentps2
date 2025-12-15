from opentps.core.processing.planOptimization.objectives.dosimetricObjectives._dosimetricObjective import DosimetricObjective
import numpy as np
import scipy.sparse as sp
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

class DUniform(DosimetricObjective):
    """
    The DUniform objective defines a constraint on the uniformity of dose distribution within a specified ROI. The penalty is determined by the variance of voxel doses relative to the prescribed uniform dose level, increasing with the magnitude of deviation from this level.

    The objective is met when all voxel doses within the ROI are equal to the specified uniform dose level.

    Attributes
    ----------
    roi : image3D
        The region of interest for which the DUniform objective is calculated.
    limitValue : int
        The Dose value that defines the target uniform dose level.
    weight : float, optional
        The weight of the objective function, default is 1.
    robust : bool, optional
        If True, the objective function is calculated in a robust way, default is False.
    """
    def __init__(self,roi,limitValue:int ,weight = 1,robust = False):
        super().__init__(metric='DUniform',roi=roi,weight=weight,robust=robust)
        if limitValue < 0:
            raise ValueError("DUniform limitValue must be greater or equal to 0 but is currently set to {}".format(limitValue))
        self.limitValue = limitValue

    def _eval(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")
        if self.GPU_acceleration:
            f = cp.mean((dose[self.maskVec_GPU] - self.limitValue) ** 2)
            self.fValue = f
        else:
            f = np.mean((dose[self.maskVec] - self.limitValue) ** 2)
            self.fValue = f
        return self.fValue

    def _grad(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")

        if self.GPU_acceleration:
            dfdD = 2 / cp.sum(self.maskVec_GPU) * (dose[self.maskVec_GPU] - self.limitValue)
        else:
            dfdD = 2 / np.sum(self.maskVec) * (dose[self.maskVec] - self.limitValue)

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