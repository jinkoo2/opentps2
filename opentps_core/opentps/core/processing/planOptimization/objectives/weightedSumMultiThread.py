from opentps.core.processing.planOptimization.objectives.baseFunction import BaseFunc
import numpy as np
import threading
import time


class WeightedSumMultiThread(BaseFunc):
    def __init__(self, Nthreads=10):
        super(WeightedSumMultiThread, self).__init__()
        self.functionList = []
        self.fValue = 0.0
        self.gradVector = None
        self.Nthreads = Nthreads

    def _eval(self, x, **kwargs):
        threads = []

        # Launch evaluation threads
        for obj in self.functionList:
            while threading.active_count() > self.Nthreads:
                time.sleep(0.01)
            t = threading.Thread(target=obj.eval, args=(x,), kwargs=kwargs)
            threads.append(t)
            t.start()

        # Wait for results and accumulate
        f = 0.0
        for idx, obj in enumerate(self.functionList):
            threads[idx].join()
            f += obj.weight * obj.fValue

        self.fValue = f
        return f

    def _grad(self, x, **kwargs):
        return_dfdD = kwargs.get('return_dfdD', False)

        if return_dfdD:
            grad = np.zeros(kwargs['dose'].shape, dtype=np.float32)
        else:
            grad = np.zeros(len(x), dtype=np.float32)

        # Launch gradient threads
        threads = []
        for obj in self.functionList:
            while threading.active_count() > self.Nthreads:
                time.sleep(0.01)
            t = threading.Thread(target=obj.grad, args=(x,), kwargs=kwargs)
            threads.append(t)
            t.start()

        # Wait for results and accumulate
        for idx, obj in enumerate(self.functionList):
            threads[idx].join()
            if return_dfdD:
                grad[obj.maskVec] += obj.weight * obj.gradVector
            else:
                grad += obj.weight * obj.gradVector

        self.gradVector = grad
        return grad