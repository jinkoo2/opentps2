from __future__ import annotations

__all__ = ['ObjectivesList']

import logging
from typing import Sequence, Union
from opentps.core.data.images._roiMask import ROIMask
from opentps.core.processing.planOptimization.objectives.baseFunction import BaseFunc
logger = logging.getLogger(__name__)


class ObjectivesList:
    """
    This class is used to store the objectives of a plan.
    A plan can have multiple objectives.
    An objective can be a Fidelity Objective or an Exotic Objective.

    Attributes
    ----------
    objectivesList: Sequence[BaseFunc]
        list of all objectives (Robust and non-robust)
    nonRobustObjList: Sequence[BaseFunc]
        list of non-robust objectives
    robustObjList: Sequence[BaseFunc]
        list of robust objectives
    targetName: str
        name of the target
    targetPrescription: float
        prescription dose of the target
    targetMask: ROIMask
        mask of the target
    """
    def __init__(self):
        self.objectivesList:Sequence[BaseFunc] = []
        self.nonRobustObjList:Sequence[BaseFunc] = []
        self.robustObjList:Sequence[BaseFunc] = []
        self.targetName:Union[str,Sequence[str]] = []
        self.targetPrescription:Union[float,Sequence[float]] = []
        self.targetMask:Union[ROIMask,Sequence[ROIMask]] = []

    def setTarget(self, roiName, roiMask, prescription):
        """
        Set the targets name and prescription doses (primary + secondary).

        Parameters
        ----------
        roiName: str
            name of the target 
        roiMask: ROIMask
            mask of the target 
        prescription: float
            prescription dose of the target
        """

        self.targetName.append(roiName)
        self.targetMask.append(roiMask)
        self.targetPrescription.append(prescription)


    def addObjective(self, function:BaseFunc):
        """
        Add an Exotic Objective to the list.

        Parameters
        ----------
        weight: float
            weight of the objective
        obj: ExoticObjective
            objective to append
        """
        if isinstance(function,BaseFunc):
            self.objectivesList.append(function)
            if function.robust:
                self.robustObjList.append(function)
            else:
                self.nonRobustObjList.append(function)

        else:
            raise ValueError(function.__class__.__name__ + ' is not a valid type for objective')