from opentps.core.processing.planOptimization.objectives.dosimetricObjectives._dosimetricObjective import DosimetricObjective
from opentps.core.data.images._image3D import Image3D
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

class EUDMax(DosimetricObjective):
    def __init__(self, roi:Image3D, limitValue:float, EUDa: float, weight = 1, robust = False):
        """
        The EUDmax objective defines a constraint on the EUD of a specified ROI.
         It penalizes cases where the EUD exceeds the prescribed maximum dose
         level. The penalty increases proportionally with the amount by which
         the EUD is greater than the threshold.

        The objective is met when the EUD of the ROI is less than or equal
        to the specified maximum dose.

        Attributes
        ----------
        roi : image3D
            The region of interest for which the EUDMax objective is calculated.
        limitValue : float
            The Dose value above which the EUD is penalized.
        EUDa : float
            The EUD parameter that defines the dose-response relationship.
        weight : float, optional
            The weight of the objective function, default is 1.
        robust : bool, optional
            If True, the objective function is calculated in a robust way, default is False.
        """
        super().__init__(metric='EUDMax',roi=roi,weight=weight,robust=robust)
        if limitValue < 0:
            raise ValueError("EUDMax limitValue must be greater or equal to 0 and is currently set to {}".format(limitValue))
        self.limitValue = limitValue
        self.EUDa = EUDa

    def _eval(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")
        if self.GPU_acceleration:
            DVH_a = cp.power(cp.mean(cp.power(dose[self.maskVec_GPU], self.EUDa)), (1 / self.EUDa))
            f = cp.power(cp.maximum(0, DVH_a - self.limitValue), 2)
            self.fValue = f
        else:
            DVH_a = np.pow(np.mean(np.pow(dose[self.maskVec], self.EUDa)),(1 / self.EUDa))
            f = (np.maximum(0, DVH_a - self.limitValue))**2
            self.fValue = f
        return self.fValue

    def _grad(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")

        if self.GPU_acceleration:
            dEUDadD = 1/cp.sum(self.maskVec_GPU) * cp.power(dose[self.maskVec_GPU], self.EUDa-1) * cp.power(cp.mean(cp.power(dose[self.maskVec_GPU], self.EUDa)), (1/self.EUDa) - 1)
            EUD_a = cp.power(cp.mean(cp.power(dose[self.maskVec_GPU], self.EUDa)),(1 / self.EUDa))
            dfdD = 2*dEUDadD*cp.maximum(0,EUD_a - self.limitValue)
        else:
            dEUDadD = 1/np.sum(self.maskVec) * np.power(dose[self.maskVec], self.EUDa-1) * np.power(np.mean(np.power(dose[self.maskVec], self.EUDa)), (1/self.EUDa) - 1)
            EUD_a = np.power(np.mean(np.power(dose[self.maskVec], self.EUDa)),(1 / self.EUDa))
            dfdD = 2*dEUDadD*np.maximum(0,EUD_a - self.limitValue)

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