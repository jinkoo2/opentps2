import os
import sys
import gzip
import numpy as np
import logging

from opentps.core.data._patientData import PatientData
from opentps.core.data.images._image2D import Image2D
from opentps.core.data.images._image3D import Image3D
from opentps.core.data.images._roiMask import ROIMask
from opentps.core.data.images._vectorField3D import VectorField3D

logger = logging.getLogger(__name__)

# Prefer SimpleITK for MHD when available (has useCompression flag)
try:
    import SimpleITK as sitk
    _SITK_AVAILABLE = True
except ImportError:
    sitk = None
    _SITK_AVAILABLE = False


def importImageMHD(headerFile):
    """
    Import image from a pair of MHD header and binary files.

    Parameters
    ----------
    headerFile: str
        Path of the MHD header file to be imported.

    Returns
    -------
    image: image3D object
        The function returns the imported image if successfully imported, or None in case of error.
    """
    if _SITK_AVAILABLE and os.path.isfile(headerFile):
        try:
            sitk_img = sitk.ReadImage(headerFile)
            arr = np.transpose(sitk.GetArrayFromImage(sitk_img))
            inputFolder, inputFile = os.path.split(headerFile)
            fileName, _ = os.path.splitext(inputFile)
            return Image3D(
                imageArray=arr,
                name=fileName,
                origin=sitk_img.GetOrigin(),
                spacing=sitk_img.GetSpacing()
            )
        except Exception as e:
            logger.debug("SimpleITK read failed, using custom MHD reader: %s", e)

    inputFolder, inputFile = os.path.split(headerFile)
    fileName, fileExtension = os.path.splitext(inputFile)

    metaData = readHeaderMHD(headerFile)
    binaryFile = metaData["ElementDataFile"] if os.path.isabs(metaData.get("ElementDataFile", "")) else os.path.join(inputFolder, metaData["ElementDataFile"])
    image = readBinaryMHD(binaryFile, metaData)
    return image



def _dtype_to_element_type(dtype):
    """Map numpy dtype to MetaImage ElementType so CT/int images can be saved as MET_SHORT."""
    dtype = np.dtype(dtype)
    if dtype == np.float64:
        return "MET_DOUBLE"
    if dtype == np.float32:
        return "MET_FLOAT"
    if dtype in (np.int16, np.uint16):
        return "MET_SHORT"  # 16-bit; typical for DICOM CT (Hounsfield units)
    if dtype in (np.int32, np.uint32):
        return "MET_INT"
    if dtype == np.bool_:
        return "MET_BOOL"
    # default for other types (e.g. int8, int64)
    return "MET_FLOAT"


def exportImageMHD(outputPath, image, vectorField='velocity'):
    """
    Export image in MHD format (header + binary files).
    Uses SimpleITK when available (with useCompression=True); otherwise custom writer.

    Parameters
    ----------
    outputPath: str
        Path of the MHD header file that will be generated.

    image: image3D (or sub-class) object or list of image3D (e.g. a 4DCT or 4D displacement field)
        Image to be exported.
    """
    # SimpleITK path: single Image3D (no list, no VectorField3D) with compression
    if _SITK_AVAILABLE and type(image) is not list and not isinstance(image, VectorField3D):
        try:
            arr = np.asarray(image._imageArray)
            # CT from DICOM is often loaded as float32; cast to int16 when values fit so we get MET_SHORT
            if arr.dtype in (np.float32, np.float64):
                arr_min, arr_max = np.nanmin(arr), np.nanmax(arr)
                if arr_min >= -32768 and arr_max <= 32767:
                    arr = arr.astype(np.int16)
            sitk_img = sitk.GetImageFromArray(np.transpose(arr))
            sitk_img.SetOrigin(np.array(image._origin, float))
            sitk_img.SetSpacing(np.array(image._spacing, float))
            mhd_path = outputPath if outputPath.lower().endswith(".mhd") else outputPath + ".mhd"
            sitk.WriteImage(sitk_img, mhd_path, useCompression=True)
            logger.info("Write MHD file (SimpleITK, compressed): " + mhd_path)
            return
        except Exception as e:
            logger.debug("SimpleITK write failed, using custom MHD writer: %s", e)

    # Custom path (list / VectorField3D / or SimpleITK unavailable)
    destFolder, destFile = os.path.split(outputPath)
    fileName, fileExtension = os.path.splitext(destFile)
    if fileExtension == ".mhd" or fileExtension == ".MHD":
      mhdFile = destFile
      rawFile = fileName + ".raw"
    else:
      mhdFile = destFile + ".mhd"
      rawFile = destFile + ".raw"
    mhdPath = os.path.join(destFolder, mhdFile)
    rawPath = os.path.join(destFolder, rawFile)

    metaData = generateDefaultMetaData()
    if type(image) is list:
        if isinstance(image[0], VectorField3D):
            metaData["ElementNumberOfChannels"] = image[0]._imageArray.shape[3]
        image = ListOfImagesToNPlus1DimImage(image)
    else:
        if isinstance(image, VectorField3D):
            metaData["ElementNumberOfChannels"] = image._imageArray.shape[3]

    metaData["NDims"] = len(image._spacing)
    metaData["DimSize"] = tuple(image.gridSize)
    metaData["ElementSpacing"] = tuple(image._spacing)
    metaData["Offset"] = tuple(image._origin)
    metaData["ElementDataFile"] = rawFile
    metaData["ElementType"] = _dtype_to_element_type(image._imageArray.dtype)

    binaryData = np.asarray(image._imageArray)
    # CT from DICOM is often float32; use MET_SHORT when values fit in int16
    if metaData["ElementType"] == "MET_FLOAT" and binaryData.dtype in (np.float32, np.float64):
        arr_min, arr_max = np.nanmin(binaryData), np.nanmax(binaryData)
        if arr_min >= -32768 and arr_max <= 32767:
            metaData["ElementType"] = "MET_SHORT"
            binaryData = binaryData.astype(np.int16)
    writeHeaderMHD(mhdPath, metaData=metaData)
    writeBinaryMHD(rawPath, binaryData, metaData=metaData)




def generateDefaultMetaData():
    """
    Generate a Python dictionary with default values for MHD header parameters.

    Returns
    -------
    metaData: dictionary
        The function returns a Python dictionary with default MHD header information
    """

    return {
        "ObjectType": "Image",
        "NDims": 3,
        "ElementNumberOfChannels": 1,
        "DimSize": (0,0,0),
        "ElementSpacing": (1.0, 1.0, 1.0),
        "Offset": (0.0, 0.0, 0.0),
        "BinaryData": True,
        "CompressedData": False,
        "ElementByteOrderMSB": False,
        "ElementType": "MET_FLOAT",
        "ElementDataFile": ""
    }



def readHeaderMHD(headerFile):
    """
    Read and parse the MHD header file.

    Parameters
    ----------
    headerFile: str
        Path of the MHD header file to be loaded.

    Returns
    -------
    metaData: dictionary
        The function returns a Python dictionary with the header information
    """

    # Parse file path
    folderPath, inputFile = os.path.split(headerFile)
    fileName, fileExtension = os.path.splitext(inputFile)

    metaData = None
    
    with open(headerFile, 'r') as fid:

        # default meta data
        metaData = generateDefaultMetaData()

        # parse header data
        for line in fid:
        
            # remove comments
            if line[0] == '#': continue
            line = line.split('#')[0]
          
            # clean the string and extract key & value
            line = line.replace('\r', '').replace('\n', '').replace('\t', ' ')
            line = line.split('=')
            key = line[0].replace(' ', '')
            value = line[1].split(' ')
            value = list(filter(len, value))

          
            if "ObjectType" in key:
                metaData["ObjectType"] = value[0]

            elif "NDims" in key:
                metaData["NDims"] = int(value[0])
            
            elif "ElementNumberOfChannels" in key:
                metaData["ElementNumberOfChannels"] = int(value[0])
            
            elif "DimSize" in key:
                # metaData["DimSize"] = (int(value[0]), int(value[1]), int(value[2]))
                metaData["DimSize"] = tuple([int(v) for v in value])
            
            elif "ElementSpacing" in key:
                # metaData["ElementSpacing"] = (float(value[0]), float(value[1]), float(value[2]))
                metaData["ElementSpacing"] = tuple([float(v) for v in value])
            
            elif "Offset" in key:
                # metaData["Offset"] = (float(value[0]), float(value[1]), float(value[2]))
                metaData["Offset"] = tuple([float(v) for v in value])
            
            elif "TransformMatrix" in key:
                # metaData["TransformMatrix"] = (float(value[0]), float(value[1]), float(value[2]), \
                #                                  float(value[3]), float(value[4]), float(value[5]), \
                #                                  float(value[6]), float(value[7]), float(value[8]))
                metaData["TransformMatrix"] = tuple([float(v) for v in value])
            
            elif "CenterOfRotation" in key:
                # metaData["CenterOfRotation"] = (float(value[0]), float(value[1]), float(value[2]))
                metaData["CenterOfRotation"] = tuple([float(v) for v in value])
          
            elif "BinaryData" in key:
                metaData["BinaryData"] = bool(value[0])
          
            elif "CompressedData" in key:
                metaData["CompressedData"] = bool(value[0])
          
            elif "ElementByteOrderMSB" in key:
                metaData["ElementByteOrderMSB"] = bool(value[0])
            
            elif "ElementType" in key:
                metaData["ElementType"] = value[0]
            
            elif "ElementDataFile" in key:
                if os.path.isabs(value[0]):
                    metaData["ElementDataFile"] = value[0]
                else:
                    metaData["ElementDataFile"] = os.path.join(folderPath, value[0])

    return metaData



def readBinaryMHD(inputPath, metaData=None):
    """
    Read and the MHD binary file.

    Parameters
    ----------
    inputPath: str
        Path of the input binary file

    metaData: dictionary
        Python dictionary with the MHD header information

    Returns
    -------
    image: image3D object
        The function returns the imported image if successfully imported, or None in case of error.
    """

    # Parse file path
    folderPath, outputFile = os.path.split(inputPath)
    fileName, fileExtension = os.path.splitext(outputFile)

    if not os.path.isfile(inputPath):
        logger.error("ERROR: file " + inputPath + " not found!")
        return None

    if metaData == None:
        metaData = generateDefaultMetaData()

    # Read binary (use inputPath; support gzip when CompressedData True)
    def _read_bytes(path, compressed):
        if compressed:
            with gzip.open(path, "rb") as fid:
                return fid.read()
        with open(path, "rb") as fid:
            return fid.read()

    raw_bytes = _read_bytes(inputPath, metaData.get("CompressedData", False))

    if metaData["ElementType"] == "MET_DOUBLE":
        data = np.frombuffer(raw_bytes, dtype=np.float64)
    elif metaData["ElementType"] == "MET_BOOL":
        data = np.frombuffer(raw_bytes, dtype=bool)
    elif metaData["ElementType"] == "MET_SHORT":
        data = np.frombuffer(raw_bytes, dtype=np.int16)
    elif metaData["ElementType"] == "MET_INT":
        data = np.frombuffer(raw_bytes, dtype=np.int32)
    else:
        data = np.frombuffer(raw_bytes, dtype=np.float32)


    if metaData["ElementNumberOfChannels"] == 1:
        data = data.reshape(metaData["DimSize"], order='F')
        if metaData["NDims"]==4:
            image = []
            for d in range(data.shape[3]):
                image.append(Image3D(imageArray=data[:,:,:,d], name=fileName, origin=metaData["Offset"][:-1], spacing=metaData["ElementSpacing"][:-1]))
        elif metaData["NDims"]==3:
            image = Image3D(imageArray=data, name=fileName, origin=metaData["Offset"], spacing=metaData["ElementSpacing"])
        elif metaData["NDims"]==2:
            image = Image2D(imageArray=data, name=fileName, origin=metaData["Offset"], spacing=metaData["ElementSpacing"])
        else:
            raise NotImplementedError(f'Method not implemented for image of size {metaData["DimSize"]}')
    else:
        if metaData["NDims"]==4:
            # data = data.reshape(np.append(metaData["DimSize"][:-1], [metaData["ElementNumberOfChannels"],metaData["DimSize"][-1]]), order='F')
            data = data.reshape(np.append(metaData["ElementNumberOfChannels"], metaData["DimSize"]), order='F')
            data = np.moveaxis(data, 0,3)
            image = []
            for d in range(data.shape[4]):
                image.append(VectorField3D(imageArray=data[:,:,:,:,d], name=fileName, origin=metaData["Offset"][:-1], spacing=metaData["ElementSpacing"][:-1]))
        elif metaData["NDims"]==3:
            data = data.reshape(np.append(metaData["ElementNumberOfChannels"],metaData["DimSize"]), order='F')
            data = np.moveaxis(data, 0,3)
            image = VectorField3D(imageArray=data, name=fileName, origin=metaData["Offset"], spacing=metaData["ElementSpacing"])
        else:
            raise NotImplementedError(f'Method not implemented for vector field of size {metaData["DimSize"]}')

    return image



def writeHeaderMHD(outputPath, metaData=None):
    """
    Write MHD header file.

    Parameters
    ----------
    outputPath: str
        Path of the MHD header file to be exported.

    metaData: dictionary
        Python dictionary with the header information
    """

    # Parse file path
    destFolder, destFile = os.path.split(outputPath)
    fileName, fileExtension = os.path.splitext(destFile)
    if fileExtension == ".mhd" or fileExtension == ".MHD":
      mhdFile = destFile
      rawFile = fileName + ".raw"
    else:
      mhdFile = destFile + ".mhd"
      rawFile = destFile + ".raw"
    mhdPath = os.path.join(destFolder, mhdFile)

    if metaData == None:
        metaData = generateDefaultMetaData()

    if metaData["ElementDataFile"] == "":
        metaData["ElementDataFile"] = rawFile

    # Write header file
    logger.info("Write MHD file: " + mhdPath)
    with open(mhdPath, "w") as fid:
        for key in metaData:
            fid.write(key + " = ")
            if isinstance(metaData[key], list) or isinstance(metaData[key], tuple):
                for element in metaData[key]: fid.write(str(element) + " ")
            else:
                fid.write(str(metaData[key]))
            fid.write("\n") 



def writeBinaryMHD(outputPath, data, metaData=None):
    """
    Write MHD binary file.

    Parameters
    ----------
    outputPath: str
        Path of the output binary file.

    data: Numpy array
        Numpy array with the image to be exported.

    metaData: dictionary
        Python dictionary with the MHD header information.
    """

    # Parse file path
    destFolder, destFile = os.path.split(outputPath)
    fileName, fileExtension = os.path.splitext(destFile)
    if fileExtension == ".raw" or fileExtension == ".RAW":
      rawFile = destFile
    elif fileExtension == ".gz":
      rawFile = destFile  # e.g. ct.raw.gz (when reading legacy compressed)
    else:
      rawFile = destFile + ".raw"
    rawPath = os.path.join(destFolder, rawFile)

    if metaData == None:
        metaData = generateDefaultMetaData()
        metaData["ElementDataFile"] = rawFile
      
    # convert data type to match ElementType
    if metaData["ElementType"] == "MET_DOUBLE" and data.dtype != np.float64:
      data = np.copy(data).astype(np.float64)
    elif metaData["ElementType"] == "MET_FLOAT" and data.dtype != np.float32:
      data = np.copy(data).astype(np.float32)
    elif metaData["ElementType"] == "MET_SHORT" and data.dtype not in (np.int16, np.uint16):
      data = np.copy(data).astype(np.int16)  # signed short, e.g. CT Hounsfield
    elif metaData["ElementType"] == "MET_INT" and data.dtype not in (np.int32, np.uint32):
      data = np.copy(data).astype(np.int32)
    elif metaData["ElementType"] == "MET_BOOL" and data.dtype != np.bool_:
      data = np.copy(data).astype(np.bool_)

    if metaData["ElementNumberOfChannels"]==3: # VectorField3D
        data = np.moveaxis(data, 3, 0)
    
    if data.dtype.byteorder == '>':
      data.byteswap()
    elif data.dtype.byteorder == '=' and sys.byteorder != "little":
      data.byteswap()

    # Write binary file (custom path is uncompressed; SimpleITK path uses its own compression)
    flat = data.reshape(data.size, order='F')
    raw_bytes = flat.tobytes()
    with open(rawPath, "wb") as fid:
        fid.write(raw_bytes)


class ListOfImagesToNPlus1DimImage(Image3D):
    
    def __init__(self, listOfImages:list):
        self.listOfImages = listOfImages
        self._spacing = np.append(listOfImages[0].spacing, 1)
        self._origin = np.append(listOfImages[0].origin, 0)
        self._imageArray = np.stack([image._imageArray for image in listOfImages], axis=-1)


    @property
    def gridSize(self):
        if isinstance(self.listOfImages[0], VectorField3D):
            return np.append(self.listOfImages[0].gridSize, self._imageArray.shape[-1])
        else:
            return self._imageArray.shape

    