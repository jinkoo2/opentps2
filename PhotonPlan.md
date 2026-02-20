# Photon plan structure in OpenTPS

This document lists the main plan objects and their properties for **photon** plans (IMRT/VMAT), with short explanations and examples. It is based on the OpenTPS core classes and the synthetic plan built in `_myscripts/dose_calc4.py`.

---

## 1. Hierarchy

```
PhotonPlan (RTPlan)
└── beams: list of PlanPhotonBeam
    └── beam segments: list of PlanPhotonSegment
        └── beamlets: list of Beamlet (x, y, MU)
```

- **Plan** = collection of beams (e.g. one AP beam, one LAO beam).
- **Beam** = one direction (gantry/couch/collimator), with one or more **segments** (control points).
- **Segment** = one aperture (jaws + MLC) and MU; it is converted into a grid of **beamlets** (small pencil beams) for the CCC dose engine.

---

## 2. Plan level: `PhotonPlan` (and base `RTPlan`)

| Property | Type | Description |
|----------|------|-------------|
| `name` | str | Plan label (e.g. `"3x3 open AP"`). |
| `_beams` | list | List of `PlanPhotonBeam` (use `plan.beams` to read). |
| `numberOfFractionsPlanned` | int | Number of treatment fractions (e.g. `1`). |
| `scanMode` | str | `"MODULATED"` for IMRT/VMAT (required for CCC). |
| `radiationType` | str | `"Photon"` for photon plans. |
| `modality` | str | DICOM modality (e.g. `"RT Plan IOD"`). |
| `seriesInstanceUID` | str | DICOM series UID. |
| `sopInstanceUID` | str | DICOM SOP instance UID. |
| `treatmentMachineName` | str | Linac name (optional). |
| `SAD_mm` | float | Source-to-axis distance in mm (photon; often taken from first beam). |

**Derived (read-only):**

- `numberOfBeamlets`: total beamlets over all segments.
- `beamSegments`: flat list of all segments from all beams.
- `beamlets`: flat list of all beamlets.
- `beamletMUs`: array of MU per beamlet (used by CCC).

**Example (from workspace.json):**

```json
{
  "name": "3x3 open AP",
  "numberOfFractionsPlanned": 1,
  "scanMode": "MODULATED",
  "radiationType": "Photon",
  "beams": [ ... ]
}
```

---

## 3. Beam level: `PlanPhotonBeam`

One treatment direction (one gantry/couch/collimator setup). Can have multiple segments (e.g. multiple control points in VMAT).

| Property | Type | Description |
|----------|------|-------------|
| `name` | str | Beam label (e.g. `"AP"`, `"Beam 1"`). |
| `id` | int | Beam number (1-based, for DICOM). |
| `gantryAngle_degree` | float | Gantry angle in degrees (IEC convention). |
| `couchAngle_degree` | float | Couch (table) angle in degrees. |
| `beamLimitingDeviceAngle_degree` | float | Collimator angle (often on segment). |
| `isocenterPosition_mm` | list [x,y,z] | Isocenter in mm in CT/image coordinates. |
| `SAD_mm` | float | Source-to-axis distance in mm (e.g. `1000`). |
| `xBeamletSpacing_mm` | float | Spacing of beamlets in the x direction in the beam’s-eye view (mm). |
| `yBeamletSpacing_mm` | float | Spacing of beamlets in the y direction (mm). |
| `beamType` | str | e.g. `"Static"` or `"Arc"`. |
| `scalingFactor` | float | Scale factor for MU (often `1`). |
| `seriesInstanceUID` | str | UID for the beam. |

Segments are stored in `_beamSegments`; use `beam.beamSegments` or `beam[i]` to access them.

**Example:**

```json
{
  "name": "AP",
  "id": 1,
  "gantryAngle_degree": 0.0,
  "couchAngle_degree": 0.0,
  "SAD_mm": 1000.0,
  "isocenterPosition_mm": [0.0, 0.0, 0.0],
  "xBeamletSpacing_mm": 5.0,
  "yBeamletSpacing_mm": 5.0,
  "beamType": "Static",
  "segments": [ ... ]
}
```

---

## 4. Segment level: `PlanPhotonSegment`

One control point: one aperture (jaws + MLC) and one MU value. The CCC engine converts segments into beamlets using `createBeamletsFromSegments()` (requires MLC coordinates).

| Property | Type | Description |
|----------|------|-------------|
| `x_jaw_mm` | list [x1, x2] | X jaw positions in mm (e.g. `[-15, 15]` for ±15 mm). |
| `y_jaw_mm` | list [y1, y2] | Y jaw positions in mm. |
| `Xmlc_mm` | array, shape (n_leaves, 4) | X MLC: each row is `[y_low, y_high, x_left, x_right]` in mm (beam’s-eye view). Required for beamlet creation. |
| `Ymlc_mm` | list / array | Y MLC positions (optional; CCC beamlet creation uses `Xmlc_mm`). |
| `mu` | float | Monitor units for this segment. |
| `gantryAngle_degree` | float | Gantry angle for this control point. |
| `couchAngle_degree` | float | Couch angle. |
| `beamLimitingDeviceAngle_degree` | float | Collimator angle. |
| `xBeamletSpacing_mm` | float | Beamlet spacing in x (usually same as beam). |
| `yBeamletSpacing_mm` | float | Beamlet spacing in y. |
| `isocenterPosition_mm` | list | Isocenter (usually same as beam). |
| `controlPointIndex` | int | Control point index (e.g. for VMAT). |
| `scalingFactor` | float | Segment scaling (often `1`). |

**`Xmlc_mm` format:** Each row describes one MLC leaf (or leaf pair) in the Y direction:

- `y_low`, `y_high`: extent of the leaf in Y (mm).
- `x_left`, `x_right`: extent of the opening in X (mm).

For an open field, all rows can have the same X extent (e.g. `[-15, 15]`) and contiguous Y strips covering the field.

**Example (3×3 cm open, 10 leaves):**

```json
{
  "gantryAngle_degree": 0.0,
  "couchAngle_degree": 0.0,
  "beamLimitingDeviceAngle_degree": 0.0,
  "x_jaw_mm": [-15.0, 15.0],
  "y_jaw_mm": [-15.0, 15.0],
  "mu": 100.0,
  "xBeamletSpacing_mm": 5.0,
  "yBeamletSpacing_mm": 5.0,
  "Xmlc_mm": [
    [-15.0, -12.0, -15.0, 15.0],
    [-12.0,  -9.0, -15.0, 15.0],
    [ -9.0,  -6.0, -15.0, 15.0],
    [ -6.0,  -3.0, -15.0, 15.0],
    [ -3.0,   0.0, -15.0, 15.0],
    [  0.0,   3.0, -15.0, 15.0],
    [  3.0,   6.0, -15.0, 15.0],
    [  6.0,   9.0, -15.0, 15.0],
    [  9.0,  12.0, -15.0, 15.0],
    [ 12.0,  15.0, -15.0, 15.0]
  ],
  "numberOfBeamlets": 36
}
```

Here `n_leaves = 10` and the field is 3×3 cm (±15 mm); the engine generates 36 beamlets (e.g. 6×6 at 5 mm spacing).

---

## 5. Beamlet (internal)

Each segment is discretized into small pencil beams (beamlets). A **Beamlet** has:

- `_x`, `_y`: position in the beam’s-eye plane (mm).
- `_mu`: monitor units for that beamlet.

You usually do not create beamlets by hand; they are created from the segment’s `Xmlc_mm` and jaw/MLC geometry by `createBeamletsFromSegments()` (called during `plan.simplify()` before CCC).

---

## 6. Building a plan in code (example)

Minimal example for one AP beam with one 3×3 cm open segment:

```python
from opentps.core.data.plan import PhotonPlan
from opentps.core.data.plan import PlanPhotonBeam
from opentps.core.data.plan import PlanPhotonSegment
import numpy as np

plan = PhotonPlan(name="3x3 open AP", patient=None)
plan.numberOfFractionsPlanned = 1
plan.scanMode = "MODULATED"
plan.radiationType = "Photon"

beam = PlanPhotonBeam()
beam.name = "AP"
beam.id = 1
beam.gantryAngle_degree = 0.0
beam.couchAngle_degree = 0.0
beam.SAD_mm = 1000.0
beam.isocenterPosition_mm = [0.0, 0.0, 0.0]
beam.xBeamletSpacing_mm = 5.0
beam.yBeamletSpacing_mm = 5.0

half = 15.0   # 3x3 cm
segment = beam.createBeamSegment()
segment.gantryAngle_degree = 0.0
segment.couchAngle_degree = 0.0
segment.beamLimitingDeviceAngle_degree = 0.0
segment.x_jaw_mm = [-half, half]
segment.y_jaw_mm = [-half, half]
n_leaves = 10
y_edges = np.linspace(-half, half, n_leaves + 1)
segment.Xmlc_mm = np.column_stack([
    y_edges[:-1], y_edges[1:],
    np.full(n_leaves, -half), np.full(n_leaves, half),
])
segment.mu = 100.0

plan.appendBeam(beam)
```

Then pass `plan` (with a CT) to `CCCDoseCalculator.computeDose(ct, plan)`.

---

## 7. Notes

- **Photon only:** This describes photon plans (`PhotonPlan`, `PlanPhotonBeam`, `PlanPhotonSegment`). Proton plans use different classes (`ProtonPlan`, `PlanProtonBeam`, layers/spots).
- **MLC required for CCC:** Segments must have `Xmlc_mm` set (non-empty, 2D) so that `createBeamletsFromSegments()` can run; otherwise the plan yields zero beamlets and no dose.
- **Coordinates:** All positions and extents are in **mm** in the appropriate coordinate system (e.g. isocenter in CT coordinates; jaws and MLC in beam’s-eye view).
- **workspace.json:** The script `dose_calc4.py` writes a serialized view of the plan (and related metadata) to `workspace.json` at the project root; the structure there matches the plan/beam/segment properties above (no binary data).

For full class definitions and methods, see:

- `opentps_core/opentps/core/data/plan/_rtPlan.py`
- `opentps_core/opentps/core/data/plan/_photonPlan.py`
- `opentps_core/opentps/core/data/plan/_planPhotonBeam.py`
- `opentps_core/opentps/core/data/plan/_planPhotonSegment.py`
