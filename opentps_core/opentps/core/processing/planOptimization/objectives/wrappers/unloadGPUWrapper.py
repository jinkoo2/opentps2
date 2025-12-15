from opentps.core.processing.planOptimization.objectives.baseFunction import BaseFunc
import numpy as np
try:
    import cupy as cp
    import cupyx as cpx
    cupy_available = True
except:
    cupy_available = False

class UnloadGPUWrapper(BaseFunc):
    """
    A wrapper class to offload computations from GPU to CPU.
    This class takes a function that performs computations on the GPU
    and ensures that the results are transferred back to the CPU.

    Attributes:
    --------
    func : BaseFunc
        The function to be wrapped, which performs computations on the GPU.
    """
    def __init__(self, func):
        super().__init__()
        self.func = func

    def _eval(self, x, **kwargs):
        x = cp.asarray(x)
        f = self.func.eval(x, **kwargs)
        self.fValue = f.get()
        return self.fValue

    def _grad(self, x, **kwargs):
        x = cp.asarray(x)
        grad = self.func.grad(x, **kwargs)
        grad = cp.asnumpy(grad)
        self.gradVector = grad
        return self.gradVector