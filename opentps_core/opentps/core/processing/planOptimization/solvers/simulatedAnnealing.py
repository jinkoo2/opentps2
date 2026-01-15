import numpy as np

from opentps.core.processing.planOptimization.solvers.solver import ConvexSolver

class SimulatedAnnealing(ConvexSolver):

    """
    Simulated Annealing optimizer.

    Attributes
    ----------
    T : float
        Initial temperature.
    maxiter : int
        Maximum number of iterations.
    coolingSchedule : function
        Cooling schedule function.
    acceptance : function
        Acceptance function.
    """


    def __init__(self,**kwargs):
        super(SimulatedAnnealing,self).__init__(**kwargs)
        self.params['T'] = self.params.get('T',1.0)
        self.params['maxiter'] = self.params.get('maxiter',1000)
        self.params['coolingSchedule'] = self.params.get('coolingSchedule',self.coolingSchedule)
        self.params['acceptance'] = self.params.get('acceptance',self.acceptance)
        self.params['ftol'] = 0.0
        self.current_value = 0.0


    def coolingSchedule(self,T):
        """
        Simple cooling schedule: T_new = 0.99 * T_old

        Parameters
        ----------
        T : float
            Current temperature.

        Returns
        -------
        float
            New temperature.
        """
        return T*0.99

    def acceptance(self,delta,T):
        """
        Acceptance function based on Metropolis criterion.
        Parameters
        ----------
        delta : float
            Change in objective function value.
        T : float
            Current temperature.

        Returns
        -------
        bool
            Whether to accept the new solution.

        """
        if delta < 0:
            return True
        else:
            return np.random.rand() < np.exp(-delta/T)

    def updateStrategy(self,x):
        """
        Simple update strategy: randomly perturb one element of x.
        Parameters
        ----------
        x : np.ndarray
            Current solution.
        Returns
        -------
        np.ndarray
            New solution.

        """
        i = np.random.randint(0,len(x))
        x[i] += np.random.randn()
        return x

    def _pre(self,function,x0):
        for f in function:
            self.smoothFuns.append(f)
        self.sol = x0
        self.current_value = np.sum([fun.eval(self.sol) for fun in self.smoothFuns])

    def _algo(self):
        x = self.updateStrategy(self.sol)
        ftot = np.sum([fun.eval(x) for fun in self.smoothFuns])
        delta = ftot - self.current_value
        if self.params['acceptance'](delta,self.params['T']):
            self.sol = x
        self.params['T'] = self.params['coolingSchedule'](self.params['T'])

    def _post(self):
        pass

