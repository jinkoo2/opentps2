from opentps.core.processing.planOptimization.objectives.baseFunction import BaseFunc
import logging
import numpy as np
import scipy.sparse as sp
from copy import copy

try:
    import cupy as cp
    import cupyx as cpx
    cupy_available = True
except:
    cupy_available = False

try:
    import sparse_dot_mkl
    MKL_available = True
except:
    MKL_available = False

logger = logging.getLogger(__name__)


class RobustWorstCase(BaseFunc):
    """
    The RobustWorstCase class is designed to handle robust optimization scenarios
    by evaluating the worst-case outcome across multiple scenarios. It extends the BaseFunc
    class and incorporates both robust and non-robust components into the objective function.

    Attributes:
        worstCaseIndex : int
            Index of the scenario that yields the worst-case outcome.
        nominalIndex : int (default: 0)
            Index of the nominal scenario.
        robustFunctions : list of BaseFunc
            List of robust function instances, each will be evaluated for each scenario.
        nonRobustFunction : BaseFunc or None
            Non-robust function instance, evaluated only once on the nominal scenario.
        robustfValues : np.ndarray
            Array to store the function values for each scenario.
        nonRobustfValue : float
            Function value of the non-robust component.
        fValue : float
            Overall function value, representing the worst-case scenario.
        gradVector : np.ndarray
            Gradient vector corresponding to the worst-case scenario.
        MKL_acceleration : bool (default: False)
            Flag to indicate if MKL acceleration is enabled.
        GPU_acceleration : bool (default: False)
            Flag to indicate if GPU acceleration is enabled.
        nScenarios : int
            Number of scenarios to consider in the robust optimization including the nominal one.
        savedWC : any
            Placeholder for saving the worst-case scenario details if needed.
        """

    def __init__(self,nScenarios,GPU_acceleration=False):
        super(RobustWorstCase, self).__init__(GPU_acceleration=GPU_acceleration)
        self.worstCaseIndex = 0
        self.nominalIndex = 0
        self.robustFunctions = []
        self.robustfValues = None
        self.nonRobustfValue = 0.0
        self.fValue = 0.0
        self.gradVector = None
        self.nonRobustFunction = None
        self.nScenarios = nScenarios
        self.savedWC = None
        if self.GPU_acceleration:
            self.robustfValues = cp.zeros(self.nScenarios, dtype=cp.float32)
        else:
            self.robustfValues = np.zeros(self.nScenarios, dtype=np.float32)

    def _eval(self, x, **kwargs):
        # Nominal scenario loop
        if self.nonRobustFunction is not None:
            self.nonRobustfValue = self.nonRobustFunction.eval(x)
        else:
            self.nonRobustfValue = 0.0

        for scenarioIndex in range(self.nScenarios):
            self.robustfValues[scenarioIndex] = self.robustFunctions[scenarioIndex].eval(x)
        if self.GPU_acceleration:
            self.worstCaseIndex = int(cp.argmax(self.robustfValues))
        else:
            self.worstCaseIndex = np.argmax(self.robustfValues)
        self.robustfValues[:] += self.nonRobustfValue
        self.fValue = self.robustfValues[self.worstCaseIndex]
        return self.fValue

    def _grad(self, x, **kwargs):
        if self.nonRobustFunction is not None:
            nonRobustGrad = self.nonRobustFunction.grad(x)
            robustGrad = self.robustFunctions[self.worstCaseIndex].grad(x)
            grad = nonRobustGrad + robustGrad
        else:
            grad = self.robustFunctions[self.worstCaseIndex].grad(x)
        self.gradVector = grad
        return self.gradVector


