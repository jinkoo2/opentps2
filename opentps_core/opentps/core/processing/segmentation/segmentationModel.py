from scipy.ndimage import binary_dilation, binary_erosion
from scipy.spatial import cKDTree
from scipy import special
from scipy.optimize import fsolve
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
import numpy as np
import scipy.stats as stats
from opentps.core.data.images._image3D import Image3D
import copy

class ProbabilisticModel():
    def __init__(self, distanceToConfidence_mm, confidence, extention_mm, mode = 'all') -> None:
        self.distance_mm = distanceToConfidence_mm
        self.confidence = confidence
        self.extention_mm = extention_mm
        self.mode = mode
        
    def pdf(self):
        return None
    
    def sample(self):
        return None

    def initializeProbabilityFunction(self):
        maxSample = self.extention_mm
        samples = np.linspace(0,maxSample, 100)
        probability = np.array([quad(self.pdf, x, np.inf)[0] for x in samples])
        self.probabilityInterpolated = CubicSpline(np.ravel(samples), probability)

    def probability(self,d_mm):
        prob = self.probabilityInterpolated(d_mm)
        return prob

    def __call__(self, image: Image3D):
        segmentationArray = copy.deepcopy(image.imageArray)
        spacing_mm = image.spacing
        segmentationArray = segmentationArray.astype(float)
        iterations = self.extention_mm / np.min(spacing_mm) ### Decrease for speed
        kernel = np.ones((3,3,3))
        erosion = binary_erosion(segmentationArray)
        dilation = binary_dilation(segmentationArray, structure = kernel, iterations = int(np.ceil(iterations)))
        surface = erosion - segmentationArray
        extention = dilation - segmentationArray
        surface_coordinates = np.argwhere(surface)
        kdtree = cKDTree(surface_coordinates)
        for voxel_index in np.argwhere(extention):
            d, nearest_neighbor_index = kdtree.query(voxel_index, workers = 5)
            d_mm = np.linalg.norm((voxel_index - surface_coordinates[nearest_neighbor_index]) * spacing_mm)
            segmentationArray[tuple(voxel_index)] = self.probability(d_mm)
        if self.mode == 'infiltration':
            segmentationArray -= image.imageArray
        return segmentationArray

class Gaussian(ProbabilisticModel):
    def __init__(self, mean_mm, distanceToConfidence_mm, confidence, extention_mm) -> None:
        super().__init__(distanceToConfidence_mm, confidence, extention_mm)
        self.mean_mm = mean_mm
        self.findSigma()
        self.initializeProbabilityFunction()
    
    def boundaryCondition(self, sigma):
        term1 = special.erf(self.mean_mm / (np.sqrt(2) * sigma))
        term2 = special.erf((self.distance_mm - self.mean_mm) / (np.sqrt(2) * sigma))
        return self.confidence * (term1 + 1) - term1 - term2

    def findSigma(self):
        sigma0 = (self.distance_mm - self.mean_mm) / 2
        self.sigma_mm = fsolve(self.boundaryCondition, sigma0)

    def pdf(self,x):
        return np.sqrt(2/np.pi) * 1 / (self.sigma_mm * (special.erf(self.mean_mm / (np.sqrt(2) * self.sigma_mm)) + 1)) * np.exp(-0.5 * ((x-self.mean_mm) / self.sigma_mm)**2) if x>=0 else 0

    def sample(self):
        a = 0
        b = np.inf
        samples = 1
        return stats.truncnorm.rvs((a - self.mean_mm) / self.sigma_mm, (b - self.mean_mm) / self.sigma_mm, loc=self.mean_mm, scale=self.sigma_mm, size=samples)[0]
    
class Uniform(ProbabilisticModel):
    def __init__(self, distanceToConfidence_mm, confidence, extention_mm) -> None:
        super().__init__(distanceToConfidence_mm, confidence, extention_mm)
        self.d_mm = self.distance_mm / self.confidence 
        self.initializeProbabilityFunction()

    def pdf(self,x):
        return 1 / self.d_mm if x>=0 and x < self.d_mm else 0
    
    def sample(self):
        return np.random.uniform(0, self.d_mm, 1)

class MixedModel(ProbabilisticModel):
    def __init__(self, modelGTM: ProbabilisticModel, modelCTM: ProbabilisticModel, mode = 'all') -> None:
        self.extention_mm = modelGTM.extention_mm + modelCTM.extention_mm
        self.GTM = modelGTM
        self.CTM = modelCTM
        self.mode = mode
        self.initializeProbabilityFunction()

    def integrand(self,y, x):
        return self.GTM.pdf(y) * quad(self.CTM.pdf, x-y, np.inf)[0]
    
    def sample(self):
        gtmSample = self.GTM.sample()
        ctmSample = self.CTM.sample()
        if self.mode == 'infiltration':
            return ctmSample
        else:
            return gtmSample + ctmSample

    def initializeProbabilityFunction(self):
        maxSample = self.extention_mm
        samples = np.linspace(0,maxSample, 100)
        infiltrationProbability = np.array([quad(lambda y: self.integrand(y, x), 0, x)[0] for x in samples])
        if self.mode == 'infiltration':
            self.probabilityInterpolated = CubicSpline(np.ravel(samples), infiltrationProbability)
        else:
            contourErrorProb = np.array([quad(self.GTM.pdf, x, np.inf)[0] for x in samples])
            self.probabilityInterpolated = CubicSpline(np.ravel(samples), contourErrorProb + infiltrationProbability)


