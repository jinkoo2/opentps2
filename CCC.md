# Collapsed Cone Convolution (CCC) Dose Algorithm in OpenTPS

This document describes the CCC photon dose calculation in OpenTPS: architecture, algorithm, file flow, and implementation notes (including recent findings and fixes).

---

## 1. Overview

- **Algorithm:** Collapsed Cone Convolution for megavoltage photon beams (e.g. 6 MV).
- **Role in OpenTPS:** Used for photon (IMRT/VMAT) dose when the plan modality is Photon; protons use MCsquare instead.
- **Reference:** The C++ engine is described in the codebase as using a "WiscPlan"-style approach; the docstring in `CCCDoseCalculator` refers to "Collapse Cone Convolution dose calculation algorithm using WiscPlan Engine."
- **License:** The CCC engine is distributed under the GNU General Public License (see `opentps_core/opentps/core/processing/doseCalculation/photons/licenseWiscPlan.txt`).

---

## 2. Architecture: Python Wrapper + C++ Executable

The CCC implementation is **not** a separate external product. It lives inside the OpenTPS repository as two layers:

| Layer | Location | Role |
|-------|----------|------|
| **Python** | `cccDoseCalculator.py` | Prepares inputs (CT, plan, kernels, geometry), writes batch scripts, launches the executable, reads back dose. |
| **C++ executable** | `CCC_DoseEngine/` | Standalone program that reads text/binary inputs and writes per-batch dose. Built as `CCC_DoseEngine_win.exe` (Windows) or `CCC_DoseEngine` (Linux). |

During a run you will see **CCC_DoseEngine_win.exe** (or the Linux binary) as separate process(es). The Python side starts one subprocess per batch via `subprocess.Popen` (e.g. running `CCC_simulation_batch0.bat`), which invokes the executable with paths to kernel list, geometry, beam specs, and output file.

---

## 3. Directory and File Layout

All paths below are relative to the **project root**.

```
opentps_core/opentps/core/processing/doseCalculation/photons/
├── cccDoseCalculator.py      # Python driver (CCCDoseCalculator class)
├── _utils.py                # Helpers (e.g. shiftBeamlets for robustness)
├── CCC_DoseEngine/          # C++ source and build
│   ├── defs.h               # Grid, kernel, beam structs; constants (NPHI, NTHETA, etc.)
│   ├── convolution.cpp      # main(); load kernels/geometry; loop over beamlets; call convolution
│   ├── parse_func.cpp       # load_kernels(), load_geometry(), pop_beam(); file format parsing
│   ├── terma_kerma.cpp      # TERMA and collision kerma from spectrum and depth
│   ├── terma_dose_masks.cpp  # terma_mask (insideness), dose_mask (cylinder)
│   ├── raytrace.cpp         # Ray tracing / depth calculation
│   ├── calc_deff.cpp        # Effective depth (deff)
│   ├── calc_dose.cpp        # Convolution of TERMA with polyenergetic kernel; dose grid
│   ├── make_poly.cpp        # Polyenergetic kernel from monoenergetic kernels
│   ├── util.cpp             # Utilities
│   ├── Makefile.win         # Builds CCC_DoseEngine_win.exe (g++)
│   └── makefile             # Linux build
├── Kernels_6MV/             # Precomputed kernels (e.g. kernel_header.txt, fluence.bin, mu_en.bin)
├── Kernels_differentFluence/
├── LINAC/                   # CT calibration (HU → density/material) for photons
│   ├── HU_Density_Conversion.txt
│   └── HU_Material_Conversion.txt
└── licenseWiscPlan.txt      # GPL v2
```

`CCCdoseEngineIO` (write CT, plan; read dose/beamlets) lives in `opentps_core/opentps/core/io/CCCdoseEngineIO.py`.

Simulation output is written under `ProgramSettings().simulationFolder` (e.g. `openTPS_workspace/Simulations/CCC_simulation/`), with subfolders such as `Geometry/`, `Outputs/`, `execFiles/`, and `BeamSpecs/`.

---

## 4. Algorithm Summary (C++ Engine)

From the source files:

1. **Inputs (command line):**  
   The executable expects four arguments: kernel list file, geometry list file, beam-spec batch file, output path for that batch.

2. **Geometry and density:**  
   CT is loaded as a 3D density grid (from `CT_HeaderFile.txt` and `CT.bin`). HU → density is done in Python before writing; the engine uses density and effective depth.

3. **Beamlets:**  
   Each beamlet is defined by SAD, position in the beam's-eye plane (xp, yp), field size (del_xp, del_yp), and unit vectors (ip, jp, kp) for the beam coordinate system. Python writes these in `pencilBeamSpecs_batch{N}.txt` via `CCCdoseEngineIO.writePlan()`.

4. **Masks:**  
   - **dose_mask:** Voxels inside a cylinder (radius aperture + safety margin) around the central axis; only these voxels receive dose.  
   - **terma_mask:** Same grid, but values are the **fraction of each voxel inside the beam aperture** (insideness). Edge voxels are upsampled (Mus × Nus × Qus) to compute this fraction.

5. **Effective depth (calc_deff):**  
   Ray tracing from the source through the density grid to compute effective depth (e.g. radiological path) for each voxel.

6. **TERMA and kerma (terma_kerma.cpp):**  
   - TERMA (Total Energy Released per unit MAss) and collision kerma are computed per voxel inside the terma_mask.  
   - Uses the incident spectrum (fluence and energy) and attenuation (μ, μ_en) from the monoenergetic kernels.  
   - Beam hardening correction is applied (Hoban et al. 1994 PMB): ratio (kerma/terma) at depth vs. at zero depth.

7. **Polyenergetic kernel (make_poly):**  
   Monoenergetic kernels (primary, first scatter, second scatter, multiple scatter, brem_annih) are combined into a single polyenergetic kernel used for convolution.

8. **Dose (calc_dose.cpp):**  
   - Dose = convolution of TERMA grid with the polyenergetic kernel.  
   - Inhomogeneity is handled by kernel scaling; inverse-square is applied after convolution rather than to the TERMA grid.  
   - Convolution uses discrete polar/azimuthal directions (NPHI, NTHETA from `defs.h`).  
   - Doses below a fraction of the maximum are cut off (`doseCutoffThreshold` in `defs.h`).

9. **Output:**  
   Per batch, the engine writes a sparse representation of the dose (e.g. `sparseBeamletMatrix_batch{N}.bin`). Python then reads these and combines them with the plan's beamlet MUs to produce the final dose image.

---

## 5. Python Side: High-Level Flow

1. **CCCDoseCalculator.computeDose(ct, plan)**  
   - Stores CT and plan; optionally converts CT from HU to density using `ctCalibration` (e.g. LINAC or MCsquare scanner folder).  
   - Cleans and prepares the simulation directory.  
   - **writeFilesToSimuDir():**  
     - Writes CT (header + density binary) and plan (beamlet specs per batch) via `CCCdoseEngineIO.writeCT()` and `CCCdoseEngineIO.writePlan()`.  
     - `writePlan()` calls `plan.simplify()` so segments are merged and beamlets are created from MLC coordinates.  
   - **writeExecuteCCCfile():**  
     - Creates batch scripts (e.g. `CCC_simulation_batch0.bat`) that invoke the CCC executable with paths to kernel list, geometry, beam specs, and output.  
   - **\_startCCC():**  
     - Starts one subprocess per batch (e.g. `CCC_DoseEngine_win.exe`); waits for all to finish.  
   - **\_importDose():**  
     - Reads the binary dose files and beamlet MUs; assembles the 3D dose grid and returns a `DoseImage`.  
   - Multiplies dose by `plan.numberOfFractionsPlanned` and returns.

2. **Batch size**  
   - `batchSize` (default in GUI can be large, e.g. 30) controls how many beamlets go into each batch and how many processes run in parallel.  
   - Large batch size can cause high memory use (e.g. "Could not allocate memory for terma_mask"). Reducing to 1–5 is recommended on limited RAM.

3. **CT calibration**  
   - Must be set (`doseCalculator.ctCalibration = calibration`) so that HU → density is correct.  
   - Calibration is read from a scanner folder (e.g. `DoseCalculationConfig().scannerFolder` or the photon `LINAC/` folder with `HU_Density_Conversion.txt` and `HU_Material_Conversion.txt`).

---

## 6. Plan and Beamlet Requirements

- **PhotonPlan** with at least one beam and segments that have **MLC coordinates** (X MLC: `Xmlc_mm`).  
- Each segment's `createBeamletsFromSegments()` in `_planPhotonSegment.py` needs `Xmlc_mm` to be a non-empty 2D array (rows = leaves, columns = y_low, y_high, x_left, x_right in mm).  
- If a segment has no X MLC (e.g. DICOM with only MLCY or no MLC), that segment yields no beamlets. If all segments are such, the plan has zero beamlets and the engine will write empty batches.

---

## 7. Recent Findings and Fixes

These items were identified and fixed during integration and scripting (dose_calc scripts, DICOM plans, water phantom):

1. **readDicomPlan — UnboundLocalError**  
   - **Issue:** In the photon branch, `plan` was used (e.g. `plan.modality = "RT Plan IOD"`) before being assigned in some code paths; in the "unknown SOPClassUID" branch, `plan` was referenced even though it was never set.  
   - **Fix (dicomIO.py):** Create `plan = PhotonPlan(...)` at the very beginning of the photon block; in the unknown-SOPClassUID `else` branch, log with `dcmFile` and `return None` instead of using `plan`.  
   - **Fix (dataLoader.py):** Only append the result of `readDicomPlan()` to the data list when it is not `None`.

2. **createBeamletsFromSegments — TypeError and empty MLC**  
   - **Issue:** When `Xmlc_mm` was empty, the code printed a message but then indexed `self.Xmlc_mm[0, 0]`, causing TypeError (list indices must be integers or slices, not tuple). Also, `Xmlc_mm` could be a list instead of a 2D numpy array.  
   - **Fix (_planPhotonSegment.py):** If `len(self.Xmlc_mm) == 0`, return immediately after the message. Convert to numpy with `Xmlc = np.asarray(self.Xmlc_mm)` and return if not 2D or empty. Use `Xmlc` for all indexing in that method.

3. **readDose — zero beamlets**  
   - **Issue:** When a plan had no beamlets (e.g. segments without MLC X), the engine wrote batch files with 0 beamlets. `CCCdoseEngineIO.readDose()` called `sp.hstack(sparseBeamletsDose)` on an empty list, causing IndexError.  
   - **Fix (CCCdoseEngineIO.py):** In the batch loop, if `numberOfBeamletsInBatch == 0` or `len(sparseBeamletsDose) == 0`, skip the hstack and dot product. After the loop, if `totalDose is None`, set `totalDose = np.zeros(header['NbrVoxels'], dtype=np.float32)` so a zero dose grid is returned and can be reshaped to a valid DoseImage.

4. **CCC_DoseEngine_win.exe as a separate process**  
   - **Clarification:** The engine is part of the same OpenTPS repo (C++ in `CCC_DoseEngine/`), built as a standalone executable and invoked by Python. It is not a separate project or third-party binary; parallelization is achieved by running multiple instances (one per batch).

5. **Memory and batch size**  
   - If you see "Could not allocate memory for terma_mask" or similar, reduce `batchSize` (e.g. to 1–5). Killing leftover `CCC_DoseEngine_win.exe` processes may be necessary before re-running if the simulation directory was locked (e.g. PermissionError on `shutil.rmtree`).

---

## 8. Build and Execution

- **Windows:** From `opentps_core/opentps/core/processing/doseCalculation/photons/CCC_DoseEngine/`, run `nmake -f Makefile.win` (or the configured make). Requires a C++ compiler (e.g. g++). Produces `CCC_DoseEngine_win.exe`.  
- **Linux:** Use the provided `makefile` in that directory to build `CCC_DoseEngine`.  
- **Execution:** Normally you do not run the executable by hand; `CCCDoseCalculator` writes the batch scripts and launches them. The executable expects paths to kernel list, geometry list, beam-spec file, and output path (see `convolution.cpp` `main()` and `writeExecuteCCCfile()` in `cccDoseCalculator.py`).

---

## 9. CCC Tuning for Better PDD / Depth–Dose Agreement

If the calculated percent depth dose (PDD) or profiles differ from your machine’s measured data (e.g. by a few percent), you can improve agreement by commissioning the following. Tuning is optional and depends on having measured beam data (PDD, profiles, output factors) for your linac and energy.

### 9.1 CT calibration (HU → density)

- **Role:** The CCC engine uses a **density** grid (not HU). Python converts HU → density using the calibration before writing `CT.bin`. Wrong density shifts depth dose (too much/little attenuation).
- **Where:** `opentps_core/.../photons/LINAC/HU_Density_Conversion.txt` and `HU_Material_Conversion.txt`. The script uses `readScanner(scanner_folder)`; if `DoseCalculationConfig().scannerFolder` does not contain `HU_Density_Conversion.txt`, the code falls back to this photon `LINAC/` folder.
- **Tuning:** Replace or edit the table in `HU_Density_Conversion.txt` with your institution’s HU–density curve (and material table if used). For a water phantom, ensure the HU range you use (e.g. 0 or small values) maps to density 1.0 g/cm³. Use a scanner-specific calibration for patient CTs.

### 9.2 Kernel set and spectrum (fluence, μ, μ_en)

- **Role:** TERMA and kerma use the **incident spectrum**: per-energy fluence, μ (attenuation), μ_en (mass energy-absorption). These come from the kernel directory’s `fluence.bin`, `mu_*.bin`, `mu_en.bin` (see `parse_func.cpp`). The polyenergetic kernel is built from monoenergetic kernels weighted by this spectrum. Wrong spectrum or attenuation shifts PDD (build-up, slope, range).
- **Where:** The Python side uses **`Kernels_6MV`** by default in `cccDoseCalculator.createKernelFilePath()`. An alternative set is `Kernels_differentFluence/` (same layout).
- **Tuning:**
  1. **Try the other set:** Set `dose_calculator._kernelsDirOverride = 'Kernels_differentFluence'` before `computeDose()` to use the alternative kernel set. Compare PDD.
  2. **Custom kernels:** To use a fully custom kernel set (e.g. for your beam energy or machine), add a copy of `Kernels_6MV` (or `Kernels_differentFluence`) under the same `photons/` folder, rename it (e.g. `Kernels_6MV_MyLinac`), put your `kernel_header.txt`, `fluence.bin`, `mu_*.bin`, `mu_en.bin` and kernel category files there, then set `dose_calculator._kernelsDirOverride = 'Kernels_6MV_MyLinac'` (or the full path). The engine expects the same file names and layout (see `kernel_header.txt`: Nradii, Nangles, Nenergies; and the parse logic in `parse_func.cpp`).
  3. **Spectrum/fluence:** If you have a measured or simulated spectrum for your beam, you can generate new fluence (and possibly μ, μ_en) files and place them in that kernel directory. This is the most impactful for PDD shape; fluence is read as a binary float array of length Nenergies.

**Kernel file formats (how to make your own)**  
The engine reads the kernel directory via the list produced by `createKernelFilePath()` (see `parse_func.cpp`, `load_kernels()`). All binary files are **raw little-endian 32-bit floats** (C `float`, 4 bytes per value), no header. Required files and formats:

| File | Format | Description |
|------|--------|-------------|
| **kernel_header.txt** | Text: first line ignored (comment); second line = `Nradii Nangles Nenergies` (three integers). Example: `24 48 14`. | Defines grid dimensions for all .bin arrays. |
| **radii.bin** | `Nradii` floats | Radial bin boundaries (e.g. cm), length = Nradii. |
| **angles.bin** | `Nangles` floats | Angular bin boundaries (e.g. rad), length = Nangles. |
| **energies.bin** | `Nenergies` floats | Energy values (MeV) for each monoenergetic kernel, length = Nenergies. |
| **fluence.bin** | `Nenergies` floats | Relative fluence (spectrum weight) per energy. |
| **mu.bin** | `Nenergies` floats | Linear attenuation coefficient μ per energy (e.g. cm⁻¹). |
| **mu_en.bin** | `Nenergies` floats | Mass energy-absorption coefficient μ_en per energy. |
| **primary.bin** | `Nradii × Nangles × Nenergies` floats | Primary kernel. Layout: for each energy in order, then `Nradii*Nangles` values with **radius index i as fast** (i.e. `value[i + j*Nradii]` for radius i, angle j). |
| **first_scatter.bin** | Same as primary | First-scatter kernel, same layout. |
| **second_scatter.bin** | Same as primary | Second-scatter kernel. |
| **multiple_scatter.bin** | Same as primary | Multiple-scatter kernel. |
| **brem_annih.bin** | Same as primary | Bremsstrahlung/annihilation kernel. |
| **total.bin** | Same as primary | Total kernel (required by parser; engine uses the five category files above). |

Example (Python) to write a 1D per-energy file and the header:

```python
import struct
import numpy as np

Nradii, Nangles, Nenergies = 24, 48, 14
# Write kernel_header.txt
with open("kernel_header.txt", "w") as f:
    f.write("Nradii Nangles Nenergies\n")
    f.write(f"{Nradii} {Nangles} {Nenergies}\n")

# Write fluence.bin (Nenergies floats)
fluence = np.ones(Nenergies, dtype=np.float32)  # or your spectrum
with open("fluence.bin", "wb") as f:
    f.write(fluence.tobytes())

# Write mu.bin, mu_en.bin similarly (Nenergies floats each).
# energies.bin: Nenergies floats (MeV).
# radii.bin: Nradii floats; angles.bin: Nangles floats.
# Each category .bin (primary, first_scatter, ...): Nradii*Nangles*Nenergies floats,
# order: energy0_rad0_ang0, energy0_rad1_ang0, ..., energy0_rad(Nr-1)_ang(Na-1), energy1_...
arr = np.zeros((Nenergies, Nangles, Nradii), dtype=np.float32)  # [E, angle, radius]
with open("primary.bin", "wb") as f:
    f.write(arr.transpose(0, 2, 1).flatten().tobytes())  # E, radius, angle → correct order
```

Easiest way to start a custom set: copy `Kernels_6MV` (or `Kernels_differentFluence`) to a new folder, then replace only the files you need (e.g. `fluence.bin`, `mu.bin`, `mu_en.bin` for a new spectrum) and keep dimensions and layout unchanged.

### 9.3 C++ engine constants (`defs.h`)

- **Role:** These affect scatter resolution, edge handling, and low-dose cutoff.
- **Where:** `opentps_core/.../photons/CCC_DoseEngine/defs.h`.
- **Relevant defines:**
  - **`doseCutoffThreshold`** (default `0.005`): Doses below this fraction of the maximum are set to zero. Raising it can trim scatter “tails” and slightly change normalization; lowering it keeps more low dose.
  - **`NPHI`** (default `12`), **`NTHETA`** (default `6`): Number of azimuthal and polar directions for the convolution. Increasing them can improve accuracy (especially off-axis and build-up) at the cost of runtime. `NTHETA` must divide `N_KERNEL_ANGLES` (48).
  - **`RSAFE`** (default `2.0` cm): Safety margin for the dose cylinder (dose_mask). Usually leave as is unless you see geometric cutoff.
  - **`Mus`, `Nus`, `Qus`** (default `5`): Upsampling for insideness at the aperture edge. Increasing can sharpen field edges and slightly affect penumbra.
- **Tuning:** After changing `defs.h`, rebuild the CCC executable (e.g. `nmake -f Makefile.win` in `CCC_DoseEngine/`). Try a small increase in `NPHI`/`NTHETA` first (e.g. `NTHETA 8` or `12`) and compare PDD and profiles; only then consider changing cutoff or upsample factors.

### 9.4 Workflow summary

1. **Baseline:** Run your current setup (e.g. jaw-only or simple field), export the central-axis PDD (e.g. from the script’s profile along the beam direction) and compare to measured PDD.
2. **CT/density:** Confirm HU → density for your phantom (e.g. water = 0 HU → 1.0 g/cm³). Adjust `LINAC/` calibration if needed.
3. **Kernels:** Default is `Kernels_6MV`; try `_kernelsDirOverride = 'Kernels_differentFluence'` if needed. For full commissioning, add a custom kernel set with your beam’s spectrum and rebuild fluence/μ/μ_en as needed.
4. **Engine resolution:** Optionally increase `NTHETA` (and `NPHI` if needed), rebuild, and re-compare.
5. **Normalization:** If the shape is good but the absolute value is off, check whether your comparison is normalized (e.g. to dmax or to a reference point); CCC output is in Gy per plan MU as defined by the plan’s MU and fraction count.

**Configurable kernel directory:** You can switch kernel sets without editing file paths by setting `dose_calculator._kernelsDirOverride`. Use a folder name under the photons directory (e.g. `'Kernels_differentFluence'`) or an absolute path to a custom kernel directory. Example: `dose_calculator._kernelsDirOverride = 'Kernels_differentFluence'` before `computeDose()` to use the alternative set.

---

## 10. References Within the Codebase

- **Python:** `opentps.core.processing.doseCalculation.photons.cccDoseCalculator` — `CCCDoseCalculator`, `computeDose()`, `_writeFilesToSimuDir()`, `_startCCC()`, `_importDose()`.  
- **I/O:** `opentps.core.io.CCCdoseEngineIO` — `writeCT()`, `writePlan()`, `readDose()`, `readBeamlets()`, beamlet file format.  
- **C++:** `opentps_core/.../photons/CCC_DoseEngine/convolution.cpp` (main), `terma_kerma.cpp`, `calc_dose.cpp`, `terma_dose_masks.cpp`, `defs.h`.  
- **Kernels:** `Kernels_6MV/` (default), `Kernels_differentFluence/` (under the photons folder); kernel list is built in `createKernelFilePath()`.  
- **CT calibration:** `LINAC/` (HU_Density, HU_Material); `opentps.core.io.scannerReader.readScanner()`.

This file serves as the single place to understand and maintain the CCC algorithm and its integration in OpenTPS.
