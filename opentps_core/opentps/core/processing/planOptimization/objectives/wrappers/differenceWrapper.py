import numpy as np
import logging
from opentps.core.processing.planOptimization.objectives.baseFunction import BaseFunc

logger = logging.getLogger(__name__)

class DifferenceWrapper(BaseFunc):
    """
    A wrapper class to compare the evaluations and gradients of multiple functions against a basis function.
    This class takes a basis function and a list of new functions, evaluates them at the same input,
    and logs the differences in their evaluations and gradients to separate CSV files in an "Output" directory if a name is provided.

    Attributes:
    --------
    basisFunction : BaseFunc
        The basis function to compare against.
    newFunctions : list of BaseFunc
        List of new functions to be compared with the basis
    file_eval : file object (default: None)
        File object for logging evaluation differences, or None if no name is provided.
    file_grad : file object (default: None)
        File object for logging gradient differences, or None if no name is provided.
    threshold : float
        Threshold for logging differences; differences below this value are not logged.
        For the gradient, the threshold is applied to the norm of the difference vector.
    """

    def __init__(self, basisFunction, newFunctions: list, name: str=None):
        super().__init__()
        self.basisFunction = basisFunction
        self.newFunctions = newFunctions
        self.file_eval = None
        self.file_grad = None
        if name is not None:
            self.file_eval = open("Output/evaluation_diff_" + name + ".csv", "w+")
            self.file_grad = open("Output/gradient_diff_" + name + ".csv", "w+")
            title = "F0;"
            for i in range(len(newFunctions)):
                title += "  F" + str(i + 1) + ";"
            self.file_eval.write(title + "\n")
            self.file_grad.write(title + "\n")
        self.threshold = 1e-5


    def _eval(self, x, **kwargs):
        fBase = self.basisFunction.eval(x, **kwargs)
        if len(self.newFunctions) == 0:
            self.file_eval.write("\n")
            return fBase

        # Always write baseline column as 0
        if self.file_eval is not None:
            self.file_eval.write(f"{0:.6e};")

        for func in self.newFunctions:
            f = func.eval(x, **kwargs)
            diff = np.abs(f - fBase)
            if diff > self.threshold:
                logger.warning(f"Difference {diff}")
                if self.file_eval is not None:
                    self.file_eval.write(f"  {diff:.6e};")
            else:
                if self.file_eval is not None:
                    self.file_eval.write(";")  # keep column spacing
        if self.file_eval is not None:
            self.file_eval.write("\n")
        return fBase

    def _grad(self, x, **kwargs):
        gradBase = self.basisFunction.grad(x, **kwargs)
        if len(self.newFunctions) == 0:
            if self.file_grad is not None:
                self.file_grad.write("\n")
            return gradBase

        # Always write baseline column as 0
        if self.file_grad is not None:
            self.file_grad.write(f"{0:.6e};")

        for func in self.newFunctions:
            gradNew = func.grad(x, **kwargs)
            if not np.allclose(gradBase, gradNew, atol=self.threshold):
                diff_norm = np.linalg.norm(gradBase - gradNew)
                logger.warning(f"Gradient difference norm {diff_norm}")
                if self.file_grad is not None:
                    self.file_grad.write(f"  {diff_norm:.6e};")
            else:
                if self.file_grad is not None:
                    self.file_grad.write(";")  # keep column spacing
        if self.file_grad is not None:
            self.file_grad.write("\n")
        return gradBase