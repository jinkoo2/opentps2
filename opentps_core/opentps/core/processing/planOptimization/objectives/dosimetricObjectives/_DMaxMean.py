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

class DMaxMean(DosimetricObjective):
    """
    The DMaxMean objective defines a constraint on the maximum of the
    mean dose across a specified ROI. It penalizes cases where the mean
    dose of the ROI exceeds the prescribed maximum dose level.
    The penalty increases proportionally with the amount by which
    the mean dose is greater than the threshold.

    The objective is met when the mean dose of the ROI is less than or
    equal to the specified maximum dose

    Attributes
    ----------
    roi : image3D
        The region of interest for which the DMean objective is calculated.
    limitValue : float
        The Dose value below which the DMean is penalized.
    weight : float, optional
        The weight of the objective function, default is 1.
    robust : bool, optional
        If True, the objective function is calculated in a robust way, default is False.
    """
    def __init__(self,roi,limitValue,weight = 1,robust = False):
        super().__init__(weight=weight, metric='DMaxMean',roi=roi, robust=robust)
        if limitValue >= 0:
            self.limitValue = limitValue
        else:
            logger.critical("DMaxMean limitValue must be greater or equal to 0 and is currently set to {}".format(limitValue))

    def _eval(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")
        if self.GPU_acceleration:
            f = cp.maximum(0, cp.mean(dose[self.maskVec_GPU]) - self.limitValue) ** 2
            self.fValue = f
        else:
            f = np.maximum(0, np.mean(dose[self.maskVec]) - self.limitValue) ** 2
            self.fValue = f
        return self.fValue

    def _grad(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")

        if self.GPU_acceleration:
            dfdD = (2/cp.sum(self.maskVec_GPU))*cp.maximum(0,cp.mean(dose[self.maskVec_GPU]) - self.limitValue)*cp.ones_like(dose[self.maskVec_GPU])
        else:
            dfdD = (2/np.sum(self.maskVec))*np.maximum(0,np.mean(dose[self.maskVec]) - self.limitValue)*np.ones_like(dose[self.maskVec])

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