# CCC Dose Calculation — Physical Concepts and Kernel Interpretation

This document explains the **physical meaning** of the CCC (Collapsed Cone Convolution) quantities and how they relate to radiation interactions. For algorithm details, file layout, and tuning, see **CCC.md**.

---

## 1. Kernel grid: radius, angle, energy

The kernel data are discretized in **radius**, **angle**, and **energy**. In radiation terms:

| Quantity | In the kernel data | Physical meaning |
|----------|--------------------|------------------|
| **Energy** | Bin index over incident photon energies (MeV), e.g. 0.2–6.0 MeV. | **Photon energy** for which that monoenergetic kernel was computed. The engine builds a polyenergetic kernel by weighting these by your spectrum (fluence × energy × μ). |
| **Radius** | Bin boundaries in **cm** (e.g. 0.05–60 cm). | **Distance** from the reference point (where TERMA/kerma is released, or the cone vertex) to the point where dose is deposited. So radius = “how far from the interaction (or beam axis) this dose contribution is.” |
| **Angle** | Bin boundaries in **degrees** (e.g. 3.75°–180°). | **Direction** (polar angle) of that distance: “at radius *r*, in direction θ.” It is the geometric angle of the line from the reference point to the deposition point (e.g. 0° along the beam, 90° perpendicular). It is **not** the Compton scattering angle of a single interaction — the kernel already encodes the net result of many interactions (primary + scatter). |

**Summary:** For a photon of a given **energy**, the kernel answers: “How much dose is deposited at **distance** (radius) and **direction** (angle) from the interaction point?” So: **radius = distance**, **angle = direction of deposition**, **energy = photon energy** of that monoenergetic kernel.

---

## 2. TERMA and kerma

- **TERMA** (Total Energy Released per unit MAss): energy transferred from the primary photon beam to the medium per unit mass (e.g. via Compton, pair production). Computed from the incident spectrum (fluence, energy), attenuation μ, and depth. Used as the “source” for the convolution.
- **Kerma** (Kinetic Energy Released per unit MAss): here **collision kerma** — energy transferred to charged particles that is **deposited locally** (excluding radiative losses). Computed from fluence, energy, and μ_en. The engine uses the **kerma/TERMA ratio** with depth for a beam-hardening correction (Hoban et al. 1994).
- **Dose** = convolution of the TERMA grid with the **polyenergetic kernel**. Inhomogeneities are handled by effective depth and kernel scaling; inverse-square is applied after convolution.

So: spectrum (fluence, μ, μ_en) → TERMA and kerma → convolved with kernel → dose.

---

## 3. Kernel categories (primary and scatter)

The kernel is split into five categories that are combined into one polyenergetic kernel:

| Category | Physical meaning |
|----------|------------------|
| **Primary** | Dose from primary (unscattered) photons — forward-peaked. Dominates **penumbra** and in-field fall-off. |
| **First scatter** | Dose from photons that have scattered once. |
| **Second scatter** | Dose from photons that have scattered twice. |
| **Multiple scatter** | Dose from photons that have scattered many times. |
| **Brem_annih** | Dose from bremsstrahlung and annihilation photons (secondary photons from charged particles). |

**Total** = sum of the five. Scatter categories add **broad lateral spread**; primary is **tighter** in radius/angle. So: primary drives edge sharpness; scatter drives flatness and off-axis dose.

---

## 4. Spectrum (fluence, μ, μ_en)

- **Fluence** (per energy): relative number of photons in each energy bin — the **incident spectrum**. It weights the monoenergetic kernels when building the polyenergetic kernel and enters TERMA/kerma. Changing fluence (e.g. softer vs harder spectrum) changes **PDD shape** (build-up, slope) and the **effective kernel shape** (more low energy → more scatter → broader lateral profiles).
- **μ** (linear attenuation coefficient): how strongly the medium attenuates the beam (cm⁻¹). Used for TERMA and depth dependence.
- **μ_en** (mass energy-absorption coefficient): fraction of transferred energy that is absorbed locally. Used for kerma and the beam-hardening correction.

All three are **per energy**; the engine uses the same spectrum for the whole field (no built-in off-axis softening).

---

## 5. What affects lateral (off-axis) profile shape

Lateral dose profiles (e.g. perpendicular to the central axis, flatness, penumbra) are influenced by:

- **Kernel shape:** The **radial and angular distribution** of the kernel (primary + scatter .bin files). Tighter primary → sharper penumbra; more scatter → broader tails and flatter in-field.
- **Spectrum:** **Fluence** (and μ, μ_en) change the polyenergetic kernel. Softer spectrum → more lateral scatter; harder → tighter profiles.
- **C++ constants (`defs.h`):** **NPHI**, **NTHETA** (convolution directions), **doseCutoffThreshold**, **Mus/Nus/Qus** (aperture edge), **RSAFE** (dose cylinder). See CCC.md §9.3.

**Flattening filter:** The engine does **not** apply an off-axis fluence or spectrum map. Lateral fluence variation (e.g. horn, flatness) must come from the **plan** (e.g. fluence map / segment weights per beamlet). The kernel and constants then determine how that fluence becomes **dose** laterally.

---

## 6. Convolution and collapsed cones

- Dose is computed by **convolving** the TERMA grid with the **polyenergetic kernel**.
- The “collapsed cone” method uses a discrete set of **cone directions** (polar and azimuthal angles from `NPHI`, `NTHETA`). For each voxel and each direction, the kernel is applied along that ray; contributions are summed. So the **angle** in the kernel is the **direction of the cone** (where dose is deposited relative to the TERMA point), not a single scattering angle.

---

## 7. Relation to CCC.md

- **CCC.md**: Algorithm steps, file layout, kernel file formats, tuning (CT, kernel set, defs.h), workflow, code references.
- **CCC_concepts.md** (this file): Physical meaning of radius/angle/energy, TERMA/kerma, kernel categories, spectrum, and what drives lateral profiles.

Use both for a full picture: **CCC.md** for “how it’s implemented,” **CCC_concepts.md** for “what it means physically.”
