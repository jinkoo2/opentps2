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

### 5.1 Overlaying dose on CT

To view a dose (e.g. **ccc_dose.mhd**) on top of a CT (e.g. **ct_original.mhd**):

1. **Load both** via **Patient data** → **Load data**: select the folder (e.g. `_outputs/1.dose_calc_ccc_ap_static_mlc`) or the individual `.mhd` files so that both `ct_original.mhd` and `ccc_dose.mhd` are loaded. They will appear in the **Patient data** tree under the same patient if loaded from the same folder.
2. **Set the CT as the primary (background) image:** In the Patient data tree, **double‑click** the CT item (e.g. "CT - ct_original" or "Image3D - ct_original"). The viewer will show the CT as the main image.
3. **Set the dose as the overlay (secondary image):**
   - **If the dose appears as "Dose - …"** (e.g. DICOM RT Dose or dose computed in OpenTPS): **double‑click** the dose item. It will be set as the secondary image and drawn on top of the CT.
   - **If the dose appears as a generic 3D image** (e.g. "Image3D - ccc_dose" when loaded from MHD): use **drag and drop**:
     - In the **viewer toolbar** (top of the viewer area), set the dropdown to **"Drop as secondaryImage"**.
     - **Drag** the dose item from the Patient data tree and **drop** it onto the viewer. The dose will be overlaid on the CT with a color scale (e.g. jet).
4. **Adjust the overlay:** Use the viewer’s **Tools** menu (e.g. secondary image options) to change colormap, window/level, or colorbar visibility for the dose.

### 5.2 Probing data (value at point and dose profile)

#### Value at a point (CT and dose under the cursor)

To see the **CT value** and **dose value** at the current mouse position (and the 3D position in mm):

1. Display the CT as primary and, if desired, dose as secondary overlay (see §5.1).
2. In the **main viewer toolbar** (the bar above the viewer with Settings, Independent views, etc.), click **Crosshair** to turn it on.
3. **Move the mouse** over the image. The viewer shows on the image:
   - **Top-left:** Primary (CT) **Value** and **Pos** (x, y, z in mm).
   - **Bottom-left:** Secondary (dose) **Value** and **Pos** (if a dose overlay is set).

The values update continuously as you move the cursor. The crosshair marks the current point in the slice. Scrolling to another slice and moving the mouse updates the position and values for that slice.

#### Dose (and CT) profile along a line

To plot **intensity vs. distance** along a line (e.g. CT and dose profile through a point):

1. Display the CT as primary and dose as secondary (see §5.1).
2. In the **viewer’s own toolbar** (Image viewer / Graph / DVH / 3D), click **Graph** to switch to the **Profile** view (plot + small toolbar with pencil, cross, disk icons).
3. In that profile toolbar, click **New profile** (pencil icon). This enables a line widget on the image.
4. Click **Image viewer** to switch back to the image. A **line** appears on the slice; drag its endpoints to position it where you want the profile (e.g. across the target).
5. Click **Graph** again. The **Profiles** plot shows:
   - **Gray curve:** primary image (CT) along the line.
   - **White curve:** secondary image (dose) along the line (if set).
   - Horizontal axis: distance along the line (mm); vertical: intensity.
6. **Dose only (different scales):** If CT and dose are both plotted, their scales (HU vs Gy) differ so much that the dose curve can be hard to see. Use the **Plot:** dropdown in the profile toolbar: choose **Secondary (dose) only** to plot only the dose along the line; the Y axis will scale to the dose range. Use **Primary (CT) only** for CT only, or **Show: Both** to show both again.
7. Use **Stop** (cross) in the profile toolbar to disable the line widget, and **Save** (disk) to export the profile plot as an image.

You can repeat **New profile** to add more lines; the first two plot slots are used for primary and secondary image, then contours if drawn.

### 5.3 Displaying beam geometry

To show the **beam geometry** (isocenter and beam direction) for a loaded DICOM plan (photon or proton):

1. **Load** the plan together with a CT (e.g. load a DICOM folder that contains CT, ROIs, plan, and dose).
2. **Set the CT as the primary image:** In the Patient data tree, **double‑click** the CT so the viewer shows the CT.
3. **Set the plan for display:** In the Patient data tree, **double‑click** the plan (e.g. "plan - …"). The viewer will display the beam geometry on top of the CT.
4. **What you see:**
   - **2D slice view (Image viewer):** For each beam, a **magenta sphere** at the isocenter and a **magenta line** in the beam direction (from the source toward the isocenter). Scroll through slices to see the isocenter and the line in different orientations.
   - **3D view:** Click **3D Image viewer** in the viewer toolbar to see a 3D representation of the beams (e.g. nozzle/orientation).
5. **Hide beams:** Use the viewer **Tools** menu → **RT plan** → **Reset** to clear the plan from the display.

Beam geometry is supported for both **photon** (e.g. DICOM IMRT/VMAT) and **proton** plans.

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

## 9. Dose comparison

Use this to compare two dose distributions (e.g. **ccc_dose.mhd** vs **plan_dose.mhd**): difference, absolute difference, or gamma (if implemented). The comparison is shown as a secondary image on the CT and in the DVH view.

### Steps to compare ccc_dose.mhd and plan_dose.mhd

1. **Load the doses** (and a CT if you want to view the comparison on the CT):
   - **Patient data** → **Load data** → select the folder that contains `ccc_dose.mhd`, `plan_dose.mhd`, and optionally a CT (e.g. `ct_original.mhd`). All will appear under the same patient if loaded from the same folder.
2. **Open the Dose comparison panel** in the left toolbox and select:
   - **Dose 1:** choose e.g. `ccc_dose` (or the name of the first MHD dose).
   - **Dose 2:** choose e.g. `plan_dose` (or the second).
   The dropdowns list both **DoseImage** data (e.g. from DICOM or from dose computed in OpenTPS) and **3D images** that are not CT (e.g. MHD-loaded dose files), so your `.mhd` doses will appear.
3. Click **Update!** The viewer switches to dose-comparison mode: the **difference** (Dose 1 − Dose 2) is overlaid as the secondary image. Dose 2 is resampled onto Dose 1’s grid automatically.
4. **View the result:**
   - Ensure the **Image viewer** is active (not DVH/Graph). The **primary** image is the CT (or first image); the **secondary** is the comparison map (difference by default). Use the viewer **Tools** menu if needed to adjust window/level or colormap.
   - Switch to **DVH** in the viewer toolbar to see DVH curves for both Dose 1 and Dose 2 in the same plot.
5. **Change the comparison metric:** In the viewer, open **Tools** → **Dose comparison** → **Metrics**:
   - **Difference:** Dose 1 − Dose 2 (can be negative).
   - **Absolute difference:** |Dose 1 − Dose 2|.
   - **Gamma:** (if implemented) gamma index map.
   You can also change **Colormap** and **Window level/width** for the comparison overlay from that menu.

---

## 10. Scripting

- **Scripting** panel: **New scripting window** opens an embedded Python window.
- **Select new script file:** Choose a `.py` file to edit/run in the OpenTPS environment (same process, so you can use `opentps.core` and loaded data).
- Use **Run** in the scripting window to execute the script. Useful for batch loading, custom dose analysis, or exporting data programmatically.

---

## 11. Tips

- **Workspace:** The default workspace directory (e.g. for simulations and temp files) can be set in **Program settings** (if available from the main window or menu).
- **Multiple patients:** Load data from different folders; each patient is listed in Patient data. Use **Split Patients** when saving if you want one file per patient.
- **DICOM + MHD:** You can mix DICOM and MHD in the same load (e.g. load a folder that has both); the loader recognizes each type and assigns data to the correct patient/series where possible.

For algorithm and kernel details (e.g. CCC), see **CCC.md** and **CCC_concepts.md** in the project root.
