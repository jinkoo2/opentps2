# StereoPHAN Simple AP Beam Dataset

## Description

Sample dataset exported from **Varian Eclipse Treatment Planning System** containing a photon treatment plan calculated on a StereoPHAN CT phantom with a simple AP (anterior-posterior) beam configuration.

## Contents

| File Type | Count | Description |
|-----------|-------|-------------|
| CT Images | 401 slices | StereoPHAN phantom CT scan |
| RT Structure Set (RS) | 1 | Contoured structures/ROIs |
| RT Plan (RP) | 1 | Simple AP beam photon plan |
| RT Dose (RD) | 1 | Calculated dose distribution from Eclipse |

## File Details

- **CT Series**: `CT.1.2.840.113619.2.278.3.481037834.151.1722371073.735.*.dcm`
- **Structure Set**: `RS.1.2.246.352.71.4.748058666681.418721.20240801125204.dcm`
- **RT Plan**: `RP.1.2.246.352.71.5.748058666681.1295230.20260202113503.dcm`
- **RT Dose**: `RD.1.2.246.352.71.7.748058666681.2077811.20260202113625.dcm`

## Source

- **TPS**: Varian Eclipse Treatment Planning System
- **Export Date**: 2026-02-02
- **Phantom**: StereoPHAN (Stereotactic Phantom)

## Usage in OpenTPS

### Option 1: Using the GUI (Load Data Button)

1. Start OpenTPS: `.\start_opentps_anaconda_windows.bat`
2. In the **left panel**, look for the **"Patient Data"** section
3. Click the **"Load data"** button at the bottom
4. Navigate to this folder and select it (you can select multiple files or the folder)
5. OpenTPS will automatically import all DICOM data (CT, Structures, Plan, Dose)

### Option 2: Using the Import Script

Run the provided script from the `_myscripts` folder:

```python
# In OpenTPS Scripting Panel:
from opentps.core.io.dataLoader import loadData

DATA_PATH = r"C:\Users\jkim20\Desktop\projects\opentps\_mydata\eclipse_tps\stereophan_simple_ap_beam"
loadData(mainWindow.patientList, DATA_PATH)
```

### Option 3: Standalone Python

```bash
cd C:\Users\jkim20\Desktop\projects\opentps
python _myscripts\import_stereophan_data.py
```

## Notes

- This dataset is for testing and development purposes
- The dose distribution was calculated using Eclipse's AAA (Analytical Anisotropic Algorithm)
- Can be used to compare OpenTPS CCC dose calculation against Eclipse results
