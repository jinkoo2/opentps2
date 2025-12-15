from opentps.core.processing.planOptimization.objectives.baseFunction import BaseFunc
import numpy as np
import logging
logger = logging.getLogger(__name__)

try:
    import cupy as cp
    import cupyx as cpx
    cupy_available = True
except:
    cupy_available = False

class WeightedSum(BaseFunc):
    def __init__(self):
        super(WeightedSum, self).__init__()
        self.functionList = []
        self.fValue = 0.0
        self.gradVector = None

    def _eval(self, x, **kwargs):
        F = 0.0
        for i in range(len(self.functionList)):
            f_i = self.functionList[i].weight * self.functionList[i].eval(x, **kwargs)
            F += f_i
        self.fValue = F
        return self.fValue

    def _grad(self, x, **kwargs):
        if kwargs.get('return_dfdD', False):
            if self.GPU_acceleration:
                dFdD = cp.zeros(kwargs['dose'].shape, dtype=cp.float32)
            else:
                dFdD = np.zeros(kwargs['dose'].shape, dtype=np.float32)

            for i in range(len(self.functionList)):
                df_idD = self.functionList[i].weight * self.functionList[i].grad(x, **kwargs)
                if self.GPU_acceleration:
                    dFdD[self.functionList[i].maskVec_GPU] += df_idD
                else:
                    dFdD[self.functionList[i].maskVec] += df_idD
            self.gradVector = dFdD

        else:
            if self.GPU_acceleration:
                dFdx = cp.zeros(len(x), dtype=cp.float32)
            else:
                dFdx = np.zeros(len(x), dtype=np.float32)

            for i in range(len(self.functionList)):
                df_idx = self.functionList[i].weight * self.functionList[i].grad(x, **kwargs)
                dFdx += df_idx
            self.gradVector = dFdx
        return self.gradVector