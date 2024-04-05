from __future__ import annotations

__all__ = ['PlanPhotonBeam']


import copy
from typing import Optional, Sequence, Union
import unittest
import sys
import numpy as np
from itertools import compress
from typing import TYPE_CHECKING
from opentps.core.data.plan._planPhotonSegment import BeamSegment
import math

class PlanPhotonBeam:
    def __init__(self):
        self._beamSegments: Sequence[BeamSegment] = []
        self.name = ""
        self.isocenterPosition_mm = []
        self.gantryAngle_degree = 0.0
        self.couchAngle_degree = 0.0
        self.id = 0
        self.scalingFactor = 1
        self.seriesInstanceUID = ""
        self.beamType = "Static"
        self.xBeamletSpacing_mm = 1
        self.yBeamletSpacing_mm = 1
        

    def __getitem__(self, segmentNb):
        return self._beamSegments[segmentNb]

    def __len__(self):
        return len(self._beamSegments)

    def __str__(self):
        s = ''
        s+= '\t\tBeam type: '+self.beamType+'\n'
        s+='\t\tThere are {} beam segments. In total they have {} beamlets.\n'.format(len(self),self.numberOfBeamlets)
        # for segment in self._beamSegments:
        #     s += 'Beam Segment\n'
        #     s += str(segment)
        return s

    @property
    def beamSegments(self) -> Sequence[BeamSegment]:
        return [beamSegment for beamSegment in self._beamSegments]

    @property
    def numberOfBeamlets(self) -> float:
        return np.sum([len(beamSegment) for beamSegment in self._beamSegments])

    def appendBeamSegment(self, segment: BeamSegment):
        self._beamSegments.append(segment)

    def createBeamSegment(self):
        segment = BeamSegment()
        segment.isocenterPosition_mm = self.isocenterPosition_mm
        segment.couchAngle_degree = self.couchAngle_degree
        segment.gantryAngle_degree = self.gantryAngle_degree
        segment.seriesInstanceUID = self.seriesInstanceUID
        segment.xBeamletSpacing_mm = self.xBeamletSpacing_mm
        segment.yBeamletSpacing_mm = self.yBeamletSpacing_mm
        self.appendBeamSegment(segment)
        return segment

    def removeBeamSegment(self, beamSegment: Union[BeamSegment, Sequence[BeamSegment]]):
        if isinstance(beamSegment, Sequence):
            beamSegments = beamSegment
            for segment in beamSegments:
                self.removeBeamSegment(segment)
            return
        self._beamSegments.remove(beamSegment)

    @property
    def beamletMUs(self):
        mu = np.array([])
        for segment in self._beamSegments:
            mu = np.concatenate((mu, segment.beamletMUs))
        return mu
  
    @beamletMUs.setter
    def beamletMUs(self, mu: Sequence[float]):
        mu = np.array(mu)
        ind = 0
        for segment in self._beamSegments:
            segment.beamletMUs = mu[ind:ind + len(segment.beamletMUs)]
            ind += len(segment.beamletMUs)

    @property
    def beamletsXY_mm(self) -> np.ndarray:
        xy = np.array([])
        for segment in self._beamSegments:
            segmentXY = list(segment.beamletsXY_mm)
            if len(segmentXY) <= 0:
                continue
            if len(xy) <= 0:
                xy = segmentXY
            else:
                xy = np.concatenate((xy, segmentXY))
        return xy
    
    @property
    def meterset(self) -> float:
        return np.sum(np.array([segment.meterset for segment in self._beamSegments]))

    @property
    def numberOfSpots(self) -> int:
        return np.sum(np.array([len(segment) for segment in self._beamSegments]))

    # @property
    # def xBeamletSpacing_mm(self):
    #     return self._beamSegments[0].xBeamletSpacing_mm
    
    # @xBeamletSpacing_mm.setter
    def setXBeamletSpacing_mm(self, xSpacing):
        self.xBeamletSpacing_mm = xSpacing
        for segment in self._beamSegments:
            segment.xBeamletSpacing_mm = xSpacing

    # @property
    # def yBeamletSpacing_mm(self):
    #     return self._beamSegments[0].yBeamletSpacing_mm
    
    # @yBeamletSpacing_mm.setter
    def setYBeamletSpacing_mm(self, ySpacing):
        self.yBeamletSpacing_mm = ySpacing
        for segment in self._beamSegments:
            segment.yBeamletSpacing_mm = ySpacing

    def copy(self):
        return copy.deepcopy(self)

    def createEmptyBeamWithSameMetaData(self):
        beam = self.copy()
        beam._beamSegments = []
        return beam

    def createBeamletsFromSegments(self):
        for i, segment in enumerate(self._beamSegments):
            print('\tCreating beamlets for beam segment {}.'.format(i))
            segment.createBeamletsFromSegments()

    def convertInPureSegments(self):
        for segment in self._beamSegments:
            if not segment.isPure():
                "Convert Segment to pure"
                pass
    
    def isInList(self, list, key, angle, tolerance):
        if len(list) == 0:
            return np.array([False])
        dicAngleDiffs = np.array([dic[key] for dic in list]) - angle
        return np.logical_and(-tolerance < dicAngleDiffs, dicAngleDiffs <= tolerance)

    def _mergeBeamSegments(self, angleTolerance):
        list = []
        for i, segment in enumerate(self._beamSegments): ### Find mergable beam segments based on the coach and gantry angles
            if len(segment.beamlets)==0:
                segment.createBeamletsFromSegments()
            coachAngleBool = self.isInList(list, 'couchAngle', segment.couchAngle_degree, 0.01)
            gantryAngleBool = self.isInList(list, 'gantryAngle', segment.gantryAngle_degree, angleTolerance/2)
            if not coachAngleBool.any() or not gantryAngleBool.any():
                list.append({'couchAngle': segment.couchAngle_degree, 'gantryAngle': segment.gantryAngle_degree, 'indexes' : [i]})
            else:
                if not (gantryAngleBool * coachAngleBool).any():
                    list.append({'couchAngle': segment.couchAngle_degree, 'gantryAngle': segment.gantryAngle_degree, 'indexes' : [i]})
                else:
                    list[np.arange(len(list))[gantryAngleBool * coachAngleBool][0]]['indexes'].append(i)            
        
        removedElements = 0
        for element in list:
            if len(element['indexes']) <= 1:
                continue
            indexes = np.sort(element['indexes'])
            gantryAngle = 0
            for i, index in enumerate(indexes):
                if i == 0:
                    newSegment = self._beamSegments[index - removedElements]
                    newbeamletMatrix = newSegment.beamletMatrixRepresentation()
                    newSegment.removeBeamlet(newSegment.beamlets)
                    gantryAngle = self._beamSegments[index - removedElements].gantryAngle_degree
                    continue
                newbeamletMatrix += self._beamSegments[index - removedElements].beamletMatrixRepresentation() ### i-1 because we are removing segments. 
                gantryAngle += self._beamSegments[index - removedElements].gantryAngle_degree
                self.removeBeamSegment(self._beamSegments[index - removedElements])
                removedElements += 1
            newSegment.gantryAngle_degree = gantryAngle / len(indexes)
            newSegment.appendBeamlets(newbeamletMatrix)
            # self.appendBeamSegment(newSegment)

    def reorderControlPointNumber(self):
        for i,segment in enumerate(self._beamSegments):
            segment.controlPointIndex = i

    def simplify(self, threshold: float = 0.0, gantryAngleTolerance: float = 0.01):
        self._mergeBeamSegments(angleTolerance= gantryAngleTolerance)
        for segment in self._beamSegments:
            segment.simplify(threshold=threshold)
        # Remove empty layers
        self._beamSegments = [segment for segment in self._beamSegments if len(segment.beamletMUs) > 0]
        self.reorderControlPointNumber()


    


if __name__ == '__main__':
    unittest.main()
