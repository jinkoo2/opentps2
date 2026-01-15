import numpy as np
import scipy.sparse as sp
from opentps.core.processing.planOptimization.objectives.dosimetricObjectives._dosimetricObjective import \
    DosimetricObjective

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


class DMin(DosimetricObjective):
    """
    The Dmin objective defines a minimum dose constraint within a specified ROI.
    Any voxel receiving a dose lower than the prescribed minimum contributes to
    the penalty, with the penalty increasing in proportion to the amount by
    which the dose falls below the threshold.

    The objective is met when all voxel doses within the ROI are greater than
    or equal to the specified minimum dose.

    Attributes
    ----------
    roi : image3D
        The region of interest for which the DMin objective is calculated.
    limitValue : float
        The Dose value below which the DMin is penalized.
    weight : float, optional
        The weight of the objective function, default is 1.
    robust : bool, optional
        If True, the objective function is calculated in a robust way, default is False.
    """
    def __init__(self,roi,limitValue,weight = 1,robust = False):
        super().__init__(weight=weight, metric='DMin',roi=roi,robust=robust)
        if limitValue < 0:
            raise ValueError("DMin limitValue must be greater or equal to 0, but is currently set to {}".format(limitValue))
        self.limitValue = limitValue


    def _eval(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")
        if self.GPU_acceleration:
            f = cp.mean(cp.minimum(0, dose[self.maskVec_GPU] - self.limitValue) ** 2)
            self.fValue = f
        else:
            f = np.mean(np.minimum(0, dose[self.maskVec] - self.limitValue) ** 2)
            self.fValue = f
        return self.fValue

    def _grad(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")
        if self.GPU_acceleration:
            dfdD = (2/ cp.sum(self.maskVec_GPU))*cp.minimum(0, dose[self.maskVec_GPU] - self.limitValue)

        else:
            dfdD = (2 / np.sum(self.maskVec))*np.minimum(0, dose[self.maskVec] - self.limitValue)

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