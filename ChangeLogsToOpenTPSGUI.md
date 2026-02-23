# Changelog: Changes to OpenTPS GUI

Summary of modifications made to the `opentps_gui` package (outside the GUI repo’s own version control).

---

## 2026-02-02

### `opentps/gui/panels/doseComparisonPanel.py`

**Dose Comparison: allow non-DICOM doses and fix signal handler**

- **Dose-like data:** Replaced `PatientDataComboBox(DoseImage)` with standard `QComboBox` and added `_dose_like_data(patient)` so that both `DoseImage` and generic `Image3D` (e.g. doses loaded from MHD) can be selected for comparison. `_to_dose_image()` converts a selected `Image3D` to `DoseImage` via `DoseImage.fromImage3D` when needed.
- **Signal signature:** `_update_dose_combos` changed from `_update_dose_combos(self)` to `_update_dose_combos(self, _data=None)` so it can be connected to `patientDataAddedSignal` and `patientDataRemovedSignal` without `TypeError` (those signals pass a data argument).

### `opentps/gui/viewer/dataViewerComponents/profileViewer.py`

**Profile viewer: plot primary or secondary only**

- **Plot mode:** Added a "Plot:" dropdown with options "Show: Both", "Primary (CT) only", and "Secondary (dose) only" so that when CT and dose scales differ, the user can plot only one series (e.g. dose profile only) for clearer inspection.
- **Behavior:** `_drawImageProfile` and `_onShowModeChanged` update the plot according to the selected mode; primary or secondary curve is cleared when not shown.

### `opentps/gui/viewer/dataViewerComponents/imageViewerComponents/rtPlanLayer.py`

**2D beam geometry: photon plans**

- **BeamLayer.setBeam:** Now supports both `PlanProtonBeam` and `PlanPhotonBeam`. For photon beams, uses `beam.isocenterPosition_mm`, `beam.gantryAngle_degree`, and `beam.couchAngle_degree` with `beamAxisPointDicomFromAngles()` from `imageTransform3D` to compute isocenter and beam axis for display (magenta sphere and line in 2D views).

### `opentps/gui/viewer/dataViewerComponents/imageViewerComponents/rtplanLayer_3D.py`

**3D beam geometry: photon plans**

- **BeamLayer_3D.setBeam:** Now supports both `PlanProtonBeam` and `PlanPhotonBeam`. For photon beams, uses `gantryAngle_degree` and `couchAngle_degree` for the internal transform so beam direction is correct in the 3D viewer.

### Documentation

- **`OpenTPS_GUI_guide.md`:** Added/updated sections for overlaying dose on CT, probing (crosshair values and dose profiles), beam geometry display (photon and proton), Dose Comparison usage, and the "Save Data" / Export settings behavior.
