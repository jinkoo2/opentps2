# OpenTPS GUI — User Guide

This guide describes how to use the OpenTPS graphical interface: loading and saving data, viewing images, and running dose computation and plan design.

---

## 1. Starting the GUI

From the project root (with the appropriate environment activated):

```bash
cd opentps_gui
python main.py
```

Or from Python:

```python
import opentps.gui as opentps_gui
opentps_gui.run()
```

See `opentps_gui/README.md` for installation and dependencies (PyQt5, pydicom, numpy, vtk, etc.).

---

## 2. Main layout

- **Left:** Toolbox (collapsible panels). Click a panel name to open it.
- **Center/right:** Viewer area (2D/3D image display, dose overlay, etc.).
- **Bottom:** Status bar.

### Toolbox panels (in order)

| Panel | Purpose |
|-------|--------|
| **Patient data** | Load data, save data, select patient and data items. |
| **ROI** | Manage ROIs/contours; select which structure set to display. |
| **Plan design** | Create/edit beams and design the plan (IMRT, IMPT, etc.). |
| **Plan optimization** | Set objectives and run plan optimization. |
| **Dose computation** | Select CT, plan, ROI; run CCC (photon) or MCsquare (proton) dose. |
| **Plan evaluation** | Evaluate plans (e.g. scenarios, confidence intervals). |
| **Dose comparison** | Compare two dose distributions. |
| **Scripting** | Open Python scripting windows to run scripts in the OpenTPS environment. |

---

## 3. Loading data

1. Open the **Patient data** panel.
2. Click **Load data**.
3. In the dialog, select **files and/or folders** to load. You can:
   - Select one or more **files** (e.g. DICOM, MHD, or serialized `.p`/`.pbz2`/`.pkl`).
   - Select a **folder**; the loader will search it recursively for supported files.
4. Confirm. Data are grouped by patient and by type (CT, dose, plan, contours, etc.).

### Supported formats

- **DICOM:** CT, MRI, PET, RT dose, RT plan (photon/ion), RT structure set, registration (rigid/deformable). Place all files for a series in one folder, or select the folder.
- **MHD:** MetaImage (`.mhd` + `.raw` or `.zraw`). Select the `.mhd` file or a folder containing the pair. Compressed `.zraw` (zlib) is supported.
- **Serialized:** OpenTPS save format (`.p`, `.pbz2`, `.pkl`, `.pickle`) from **Save data**. Can contain full patient/data trees.

Loaded data appear in the **Patient data** tree. Select a patient to see their data; click an item to display it in the viewer (if applicable). Right‑click items for context actions (rename, delete, export, image info).

---

## 4. Saving data

### Save data (serialized)

1. **Patient data** panel → **Save data**.
2. Choose **folder and file name** (e.g. `Patient` → creates `Patient.p` or similar).
3. Options:
   - **Compress data:** Smaller files, slower save.
   - **Split Patients:** If you have multiple patients, save each to a separate file.
4. Confirm. This writes the current patient’s data tree in OpenTPS serialized format so you can reload it later with **Load data**.

### Export (DICOM / MHD / other)

1. **Patient data** panel → right‑click a data item → use the **export/save** option if available, **or** use the global export window (if exposed in the menu).
2. Alternatively, **Save data** opens a dialog that can lead to the **Export settings** window.
3. In **Export settings** you choose a folder, then for each data type (Image, Dose, Plan, Contours, Other) you select format: **DICOM**, **MHD**, **MCSQUARE**, or **PICKLE** (serialized). Click **Select folder and export** to write files.

So: **Save data** = quick save of the whole patient as serialized; **Export** = choose folder and per‑type format (DICOM/MHD/etc.) for the current patient.

---

## 5. Viewing data

- **Patient data:** Select a patient, then click a CT, dose, contour set, or other image-like data to display it in the viewer.
- **Viewer:** The central area shows 2D slices (and 3D where available). Use the viewer toolbar and mouse to scroll slices, zoom, adjust window/level, and toggle overlays (e.g. dose on CT, contours).
- **ROI panel:** Choose which structure set/ROIs are shown and how (e.g. which contours are visible).

---

## 6. Dose computation (photon CCC)

1. **Patient data:** Load a CT and a photon plan (and optionally contours).
2. **Dose computation** panel:
   - **CT:** Select the CT to use.
   - **Plan:** Select the plan.
   - **Overwrite outside this ROI:** Optionally choose an ROI; voxels outside it can be set to air for the calculation.
   - **Batch size:** Number of parallel CCC batches (e.g. 30).
3. Click **Compute dose**. The CCC engine runs; when finished, the calculated dose is added to the patient and can be displayed in the viewer.

Photon dose uses the CCC (Collapsed Cone Convolution) engine and the kernel set configured in the code (e.g. `Kernels_6MV` or `Kernels_My6MV`). Proton dose uses MCsquare when the selected plan is proton.

---

## 7. Plan design

1. **Plan design** panel: Create or select a plan design (IMRT or IMPT).
2. **Add beam:** Define beams (gantry, couch, isocenter, MLC/jaws, etc.) via the beam dialog.
3. **Design plan:** Run the design step to generate the deliverable plan from the beams and modality.

Beams and segments are defined here; then **Dose computation** uses that plan to compute dose.

---

## 8. Plan optimization

1. **Plan optimization** panel: Select CT, plan, and (for IMPT) optionally set robustness settings.
2. **Open objectives panel:** Define dose–volume or other objectives per ROI.
3. **Optimize plan:** Run optimization. For protons, beamlet computation and optimization are triggered from this panel; for photons, optimization uses the CCC dose engine as configured.

---

## 9. Scripting

- **Scripting** panel: **New scripting window** opens an embedded Python window.
- **Select new script file:** Choose a `.py` file to edit/run in the OpenTPS environment (same process, so you can use `opentps.core` and loaded data).
- Use **Run** in the scripting window to execute the script. Useful for batch loading, custom dose analysis, or exporting data programmatically.

---

## 10. Tips

- **Workspace:** The default workspace directory (e.g. for simulations and temp files) can be set in **Program settings** (if available from the main window or menu).
- **Multiple patients:** Load data from different folders; each patient is listed in Patient data. Use **Split Patients** when saving if you want one file per patient.
- **DICOM + MHD:** You can mix DICOM and MHD in the same load (e.g. load a folder that has both); the loader recognizes each type and assigns data to the correct patient/series where possible.

For algorithm and kernel details (e.g. CCC), see **CCC.md** and **CCC_concepts.md** in the project root.
