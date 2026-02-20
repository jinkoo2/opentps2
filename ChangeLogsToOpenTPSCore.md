# Changelog: Changes to OpenTPS Core

Summary of modifications made to the `opentps_core` package (outside the core repo’s own version control).

---

## 2026-02-20

### `opentps/core/data/plan/_planPhotonSegment.py`

**`createBeamletsFromSegments()` – handle missing jaw positions**

- **Issue:** For some IMRT DICOM plans, not every control point includes ASYMX/ASYMY jaw positions. Segments then had `x_jaw_mm` or `y_jaw_mm` empty or with fewer than two elements, causing `IndexError` when creating beamlets (e.g. `self.y_jaw_mm[0]`, `self.y_jaw_mm[1]`).
- **Change:** After validating `Xmlc_mm`, jaw bounds are inferred from the MLC when not set by DICOM:
  - If `len(self.y_jaw_mm) < 2`: set `y_jaw_mm = [min(Xmlc[:, 0]), max(Xmlc[:, 1])]` (y extent from MLC).
  - If `len(self.x_jaw_mm) < 2`: set `x_jaw_mm = [min(Xmlc[:, 2]), max(Xmlc[:, 3])]` (x extent from MLC).
- **Effect:** IMRT plans (e.g. `2.dose_calc_ccc_ap_imrt.py`) can run CCC without crashing when some segments lack explicit jaw positions.

### `opentps/core/io/CCCdoseEngineIO.py`

**`writePlan()` – segment index vs beam index for VMAT**

- **Issue:** The index from `argwhere(beamNumber < bemaletNumberPerSegmentAccumulated)` is a **segment** index (0 to num_segments-1). It was used to index `plan.beams`, which has one entry per **beam**. For VMAT (e.g. 1 beam, 65 segments), this caused `IndexError: list index out of range` when segment index ≥ 1.
- **Change:** Build `segmentIndexToBeamIndex` mapping (segment index → beam index). Use **segment index** for `planBeamSegments` and `rotatedVersorsPerBeam`; use **beam index** (from the mapping) only for `plan.beams[beamIndex]` (SAD).
- **Effect:** VMAT plans (e.g. `3.dose_calc_ccc_ant_vmat.py`) can run CCC without crashing in `writePlan`.
