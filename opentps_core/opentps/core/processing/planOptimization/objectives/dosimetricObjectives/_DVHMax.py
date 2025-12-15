from opentps.core.processing.planOptimization.objectives.dosimetricObjectives._dosimetricObjective import DosimetricObjective
import time
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

class DVHMax(DosimetricObjective):
    def __init__(self,roi,limitValue:float,volume:float ,weight = 1,robust = False):
        """
        The DVHMax objective defines a constraint within a specified ROI.
        It penalizes cases where the dose delivered to a given percentage
        of the ROI volume exceeds the prescribed maximum dose level.
        The penalty increases proportionally with the amount by which the
        DVH curve lies above the constraint line.

        The objective is met when the DVH of the ROI does not exceed
        the specified maximum dose at the defined volume fraction.

        Attributes
        ----------
        roi : Image3D
            The region of interest for which the DVHMax objective is calculated.
        limitValue : float
            The maximum dose value that the DVH should not exceed.
        volume : float
            The volume fraction of the ROI that is considered for the DVHMax calculation.

        """
        super().__init__(metric='DVHMax',roi=roi,weight=weight,robust=robust)
        if limitValue < 0:
            raise ValueError("DVHMax limitValue must be greater or equal to 0 but is set to {}".format(limitValue))
        if volume <= 0 or volume > 1:
            raise ValueError("DVHMax volume must be in the range (0, 1] but is set to {}".format(volume))
        self.limitValue = limitValue
        self.volume = volume

        self.active_region = None
        self.DaV = 0.0

    def _eval(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")
        if self.GPU_acceleration:
            self.DaV = self._calcInverseDVH(self.volume, dose[self.maskVec_GPU])
            self.active_region = (dose[self.maskVec_GPU] >= self.limitValue) & (dose[self.maskVec_GPU] <= self.DaV)
            f = cp.sum((dose[self.maskVec_GPU][self.active_region] - self.limitValue) ** 2) / cp.sum(self.maskVec_GPU)
            self.fValue = f
        else:
            self.DaV = self._calcInverseDVH(self.volume, dose[self.maskVec])
            self.active_region = (dose[self.maskVec] >= self.limitValue) & (dose[self.maskVec] <= self.DaV)
            f = np.sum((dose[self.maskVec][self.active_region] - self.limitValue) ** 2) / np.sum(self.maskVec)
            self.fValue = f
        return self.fValue


    def _grad(self, x, **kwargs):
        dose = kwargs.get('dose', None)
        if dose is None:
            logger.error("Dose must be provided")
            raise ValueError("Dose must be provided")

        if self.GPU_acceleration:
            dfdD = cp.zeros_like(dose[self.maskVec_GPU])
            dfdD[self.active_region] = 2*(dose[self.maskVec_GPU][self.active_region] - self.limitValue) /cp.sum(self.maskVec)

        else:
            dfdD = np.zeros_like(dose[self.maskVec])
            dfdD[self.active_region] = 2*(dose[self.maskVec][self.active_region] - self.limitValue) /np.sum(self.maskVec)

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