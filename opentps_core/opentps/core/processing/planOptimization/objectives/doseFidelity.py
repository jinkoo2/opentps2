import numpy as np
import scipy.sparse as sp
import logging
logger = logging.getLogger(__name__)

from opentps.core.processing.planOptimization.objectives.baseFunction import BaseFunc

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


class DoseFidelity(BaseFunc):
    """
    Dose fidelity objective class. Inherits from objectiveFidelity.

    Attributes
    ----------
    list : list
        List of objectives
    xSquared : bool
        If true, the weights are squared. If false, the weights are not squared.
    beamlets : sparse matrix
        Beamlet matrix
    scenariosBL : list
        List of scenarios
    GPU_acceleration : bool (default: False)
        If true, the GPU is used for the computation of the fidelity function and gradient.
    MKL_acceleration : bool (default: False)
        If true, the MKL is used for the computation of the fidelity function and gradient.
    dose : array or cupy array if GPU_acceleration is True
        Dose distribution
    function : objective function
        Objective function
    unionMaskVec : array or cupy array if GPU_acceleration is True
        Union mask vector for cropped multiplication
    croppedMultiplication : bool
        If true, cropped multiplication is used for the computation of the fidelity gradient.
    """
    def __init__(self, beamlets, xSquared=True, GPU_acceleration=False, MKL_acceleration=False):
        super(DoseFidelity, self).__init__(xSquared=xSquared, GPU_acceleration=GPU_acceleration, MKL_acceleration=MKL_acceleration)
        self.beamlets = beamlets
        self.dose = None
        self.function = None
        self.unionMaskVec = None
        self.croppedMultiplication = False

        if self.GPU_acceleration:
            self.beamlets = cpx.scipy.sparse.csc_matrix(beamlets)
            if self.croppedMultiplication:
                self.unionMaskVec = cp.asarray(self.unionMaskVec)


    def update_dose(self,x,beamlets):
        """
        Updates the dose distribution based on the current weights and beamlet matrix.

        Parameters
        ----------
        x : array or cupy array if GPU_acceleration is True
            Weights

        beamlets : scipy sparse matrix or cupy sparse matrix if GPU_acceleration is True
            Beamlet matrix

        Returns
        -------
        dose : array or cupy array if GPU_acceleration is True
            Updated dose distribution

        """
        if self.xSquared:
            # w = x^2
            if self.GPU_acceleration:
                w = cp.square(x.astype(cp.float32))
            else:
                w = np.square(x).astype(np.float32)
        else:
            # w = x
            if self.GPU_acceleration:
                w = x.astype(cp.float32)
            else:
                w = x.astype(np.float32)

        if self.GPU_acceleration:
            dose = cp.sparse.csc_matrix.dot(beamlets, w)

        elif self.MKL_acceleration:
            dose = sparse_dot_mkl.dot_product_mkl(beamlets, w)
        else:
            dose = sp.csc_matrix.dot(beamlets, w)
        return dose

    def dDdx(self,x,beamlets):
        """
        Computes the derivative of the dose with respect to the weights.
        dDdx = beamlets if weights are not squared
        dDdx = beamlets * dwdx if weights are squared
        where dwdx = 2*x
        Parameters
        ----------
        x : array or cupy array if GPU_acceleration is True
            Weights
        beamlets : scipy sparse matrix or cupy sparse matrix if GPU_acceleration is True
            Beamlet matrix

        Returns
        -------
        dDdx : scipy sparse matrix or cupy sparse matrix if GPU_acceleration is True
            Derivative of the dose with respect to the weights

        """
        # If the weights are squared, the change of variable is dDdx = dDdw * dwdx
        # pointwise multiplication is much faster than creating a diagonal matrix with dwdx and then multiplying
        if self.GPU_acceleration:
            dwdx = 2 * x
            dDdx = cpx.scipy.sparse.csc_matrix.multiply(beamlets, dwdx)
        else:
            dwdx = 2 * x
            dDdx = sp.csc_matrix.multiply(beamlets, dwdx).tocsc()
        return dDdx

    def computeFidelityFunction(self, x, dose):
        """
        Computes the fidelity function.

        Parameters
        ----------
        x : array
            Weights
        returnWorstCase : bool
            If true, the worst case scenario is returned. If false, the nominal scenario is returned.

        Returns
        -------
        fTot : float
            Fidelity function value
        worstCase : int (only if robust objectives are present)
            Worst case scenario index (-1 for nominal)
        """
        F = self.function.eval(x,dose=dose)
        return F


    def computeFidelityGradient(self, x,dose,dDdx):
        """
        Computes the fidelity gradient.

        Parameters
        ----------
        x : array
            Weights

        Returns
        -------
        dfTot : array
            Fidelity gradient

        Raises
        ------
        Exception
            If the objective metric is not supported.
        """
        if self.croppedMultiplication:
            dFdD = self.function.grad(x, dose=dose, dDdx=dDdx, return_dfdD=True)
            if self.GPU_acceleration:
                dFdx = cp.sparse.csc_matrix.dot(dDdx[:, self.unionMaskVec], dFdD[self.unionMaskVec])
            elif self.MKL_acceleration:
                dFdx = sparse_dot_mkl.dot_product_mkl(dDdx[:, self.unionMaskVec], dFdD[self.unionMaskVec])
            else:
                dFdx = sp.csc_matrix.dot(dDdx[:, self.unionMaskVec], dFdD[self.unionMaskVec])

        else:
            dFdx = self.function.grad(x, dose=dose, dDdx=dDdx, return_dfdD=False)
        return dFdx


    def _eval(self, x, **kwargs):
        self.dose = self.update_dose(x,self.beamlets)
        f = self.computeFidelityFunction(x,self.dose)
        self.fValue = f
        return self.fValue

    def _grad(self, x, **kwargs):
        dose = self.dose
        if self.xSquared:
            dDdx = self.dDdx(x,self.beamlets)
        else:
            dDdx = self.beamlets
        dDdx = dDdx.T
        g = self.computeFidelityGradient(x,dose,dDdx)
        self.gradVector = g
        return self.gradVector