import time
import logging
from opentps.core.processing.planOptimization.objectives.baseFunction import BaseFunc

logger = logging.getLogger(__name__)

class TimingWrapper(BaseFunc):
    """
    A wrapper class to measure and log the time taken for function evaluations and gradient computations.
    This class takes a function and records the time taken for each evaluation and gradient computation,
    writing the results to separate CSV files in an "Output" directory.

    Attributes:
    --------
    func : BaseFunc
        The function to be wrapped, which performs the actual computations.
    file_eval : file object
        File object for logging evaluation times.
    file_grad : file object
        File object for logging gradient computation times.
    """
    def __init__(self, func, name: str):
        super().__init__()
        self.func = func
        self.file_eval = open("Output/evaluation_time_" + name + ".csv", "w+")
        self.file_grad = open("Output/gradient_time_" + name + ".csv", "w+")

        # Header row (only one column now)
        self.file_eval.write("F0;\n")
        self.file_grad.write("F0;\n")

    def _eval(self, x, **kwargs):
        t0 = time.time()
        f = self.func.eval(x, **kwargs)
        t = time.time() - t0
        self.file_eval.write(f"{t:.6e};\n")
        return f

    def _grad(self, x, **kwargs):
        t0 = time.time()
        grad = self.func.grad(x, **kwargs)
        t = time.time() - t0
        self.file_grad.write(f"{t:.6e};\n")
        return grad