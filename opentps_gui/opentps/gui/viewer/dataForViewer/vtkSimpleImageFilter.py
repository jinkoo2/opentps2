from vtkmodules.util.vtkAlgorithm import vtkPythonAlgorithm
from vtkmodules.vtkCommonDataModel import vtkImageData

class VtkSimpleImageFilter(vtkPythonAlgorithm):
    def __init__(self):
        super().__init__()
        self.SetNumberOfInputPorts(1)
        self.SetNumberOfOutputPorts(1)

    def RequestData(self, request, inInfoVec, outInfoVec):
        # Get input and output vtkImageData objects
        input_image = vtkImageData.GetData(inInfoVec[0].GetInformationObject(0))
        output_image = vtkImageData.GetData(outInfoVec.GetInformationObject(0))

        # Shallow copy input to output
        output_image.ShallowCopy(input_image)
        return 1

