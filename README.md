# Java Coastline Change Analysis 1990–2025
### Automated Satellite-Derived Shoreline Extraction with Tidal Correction and DSAS-Equivalent Change Analysis

---

## Overview

This repository implements a three-stage automated framework for coastline change analysis along Java's coast from 1990 to 2025. The workflow integrates multi-sensor Landsat and Sentinel-2 imagery processed in Google Earth Engine, FES2022-based tidal and sea-level rise correction, and a DSAS-equivalent transect analysis producing EPR, LRR, and NSM metrics — all without manual intervention.

The shoreline extraction methodology is built upon the **CoastSat framework** (Vos et al., 2019), which established the foundational approach of using MNDWI-based water index thresholding, sub-pixel resolution shoreline mapping, and tidal correction from satellite imagery. This repository extends CoastSat's core principles by:

- Reimplementing the extraction pipeline natively in **Google Earth Engine JavaScript** for scalable batch processing across 8 epochs and multiple ROIs
- Replacing CoastSat's FES2014 tidal correction with **FES2022** constituents for improved accuracy in Indonesian waters
- Adding **land subsidence** and **sea-level rise** components to the vertical offset correction
- Integrating a **DSAS-equivalent** transect analysis (EPR, LRR, NSM) as a downstream processing stage

> **Reference framework:**
> Vos, K., Splinter, K.D., Harley, M.D., Simmons, J.A., Turner, I.L. (2019).
> CoastSat: A Google Earth Engine-enabled Python toolkit to extract shorelines from publicly available satellite imagery.
> *Environmental Modelling & Software*, 122, 104528.
> https://doi.org/10.1016/j.envsoft.2019.104528
>
> GitHub: https://github.com/kvos/CoastSat

This framework is adapted from methods developed in the author's Master's thesis:
>
> Putra, M.Z.R. (2025). *The Development of a Tidally Corrected Automatic
> Coastline Extraction Framework for Supporting Indonesia's Maritime Baseline
> Monitoring*. Master's Thesis, Geomatic Engineering, Universitas Gadjah Mada,
> Yogyakarta, Indonesia.
>
> **Author:** Mohammad Zulfi Rahadi Putra
>
> The original thesis applied this methodology to Indonesia's 183 official
> archipelagic basepoints under PP 38/2002, using FES2022 tidal normalization
> to produce datum-consistent shoreline positions at a national scale. This
> repository re-applies the same core framework to a regional coastline change
> study along Java (1990–2025), extending the temporal coverage to 8 epochs
> and integrating DSAS-equivalent change metrics (EPR, LRR, NSM).

---

## Repository Structure

```
├── Part1_GEE/
│   └── coastline_extraction.js       # Google Earth Engine extraction script
│
├── Part2_TidalCorrection/
│   ├── tidal_slr_correction.py       # FES2022 tidal + SLR vertical offset
│   ├── fes2022_indonesia_clipped.yaml
│   └── fes2022_resolved.yaml         # auto-generated at runtime
│
├── Part3_DSAS/
│   └── dsas_analysis.py              # EPR / LRR / NSM transect analysis
│
├── shorelines/                       # GEE-exported raw shorelines (.shp/.geojson)
├── shorelines_corrected/             # Output from Part 2 (tidal-corrected)
├── output_dsas/                      # Output from Part 3 (metrics + maps)
│
├── requirements.txt
└── README.md
```

---

## Workflow

```
[Part 1 — GEE]                [Part 2 — Python]              [Part 3 — Python]
Landsat 4/5/7/8/9          →  FES2022 Tidal Correction   →   DSAS Analysis
Sentinel-2                     SLR Vertical Offset             EPR / LRR / NSM
MNDWI + Canny Edge             Subsidence Component            Transect Maps
(CoastSat-inspired)            (FES2022 / Vos et al.)          Change Profiles
Export .shp per year           Export corrected .geojson
```

---

## Relationship to CoastSat

This framework adopts and extends the following CoastSat components:

| CoastSat Component | This Framework |
|---|---|
| MNDWI water index for water/land separation | Retained — implemented natively in GEE |
| Otsu / histogram thresholding | Replaced by fixed MNDWI threshold + Canny edge detection for batch GEE processing |
| Sub-pixel shoreline mapping | Approximated via morphological smoothing + Canny sigma tuning |
| FES2014 tidal correction | Upgraded to **FES2022** with subsidence and SLR components |
| Single-site processing | Extended to **multi-ROI batch processing** across 8 epochs |
| Python + GEE API | GEE JavaScript for extraction; Python for correction and analysis |

For single-site, interactive shoreline extraction with full sub-pixel accuracy, the original CoastSat toolkit is recommended: https://github.com/kvos/CoastSat

---

## Part 1 — Coastline Extraction (Google Earth Engine)

**File:** `Part1_GEE/coastline_extraction.js`

### What it does
- Loads Landsat 4, 5, 7, 8, 9 collections (1990–2025)
- Applies QA_PIXEL cloud masking and surface reflectance scaling
- Computes MNDWI following the CoastSat spectral index approach (Vos et al., 2019)
- Applies connected-pixel cleaning and morphological smoothing
- Detects coastline edges using Canny Edge Detection
- Exports shoreline vectors (`.shp`) and MNDWI rasters per epoch to Google Drive

### Years processed
`1990, 1995, 2002, 2005, 2010, 2015, 2020, 2025`

### Sensor mapping

| Year | Sensor | GEE Collection |
|---|---|---|
| 1990, 1995 | Landsat 4/5 TM | `LANDSAT/LT04/C02/T1_L2`, `LANDSAT/LT05/C02/T1_L2` |
| 2002, 2005, 2010 | Landsat 5/7 | `LANDSAT/LT05/C02/T1_L2`, `LANDSAT/LE07/C02/T1_L2` |
| 2015, 2020 | Landsat 8 OLI | `LANDSAT/LC08/C02/T1_L2` |
| 2025 | Landsat 8/9 OLI | `LANDSAT/LC08/C02/T1_L2`, `LANDSAT/LC09/C02/T1_L2` |

### Setup
1. Open [Google Earth Engine Code Editor](https://code.earthengine.google.com)
2. Paste the contents of `coastline_extraction.js`
3. Define your ROI geometry as `ROI8` (or rename the variable at the top)
4. Run and submit export tasks to Google Drive
5. Download exported files into the `shorelines/` folder

### Naming convention expected by Part 2
```
Shoreline_1990_Line.shp
Shoreline_1995_Line.shp
...
Shoreline_2025_Line.shp
```

---

## Part 2 — Tidal & SLR Correction (Python)

**File:** `Part2_TidalCorrection/tidal_slr_correction.py`

### What it does
Extends the CoastSat tidal correction approach (Vos et al., 2019) using FES2022 instead of FES2014, with additional vertical offset components:

$$\text{LWL}_{\text{corrected}} = \text{LAT}_{\text{FES2022}} - \Delta\text{MSL}_{\text{tidal}} - \text{SLR} - \text{Subsidence}$$

Where:
- $$\Delta\text{MSL}_{\text{tidal}}$$ — epoch MSL change from FES2022 hourly sampling
- $$\text{SLR}$$ — cumulative sea-level rise at 3.7 mm/yr from reference year
- $$\text{Subsidence}$$ — cumulative land subsidence by tidal class (Abidin et al., 2011)

### FES2022 Setup
FES2022 ocean tide model files are required and are **not included** in this repository due to file size. Download from:
> [AVISO+ FES2022](https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html)

Place the constituent NetCDF files in a local directory and update `fes2022_indonesia_clipped.yaml` to point to their paths.

### Key parameters (editable in script)

| Parameter | Default | Description |
|---|---|---|
| `SLR_RATE_MM_YR` | `3.7` | Global mean SLR rate (mm/yr) |
| `SUBSIDENCE_RATE_MM_YR` | micro: 4.5, meso: 12.0, macro: 7.5 | Land subsidence by tidal class (mm/yr) |
| `REF_YEAR` | `2002` | Reference epoch for delta MSL calculation |

### Outputs
```
shorelines_corrected/
├── shoreline_tidal_corrected_{year}.geojson
├── shoreline_tidal_corrected_{year}.shp
├── tidal_correction_all_years.csv
├── lwl_correction_overview.png
└── lwl_all_years_comparison.png
```

---

## Part 3 — DSAS-Equivalent Change Analysis (Python)

**File:** `Part3_DSAS/dsas_analysis.py`

### What it does
Implements the fundamental transect-based principles of the Digital Shoreline Analysis System (DSAS; Himmelstoss et al., 2021), extended with automated batch-processing across all epochs without manual intervention:

1. Builds a reference baseline from the earliest available shoreline
2. Generates shore-perpendicular transects at configurable spacing
3. Intersects each transect with all epoch shorelines
4. Computes three standard DSAS metrics per transect:

| Metric | Formula | Description |
|---|---|---|
| **NSM** | $$d_{\text{last}} - d_{\text{first}}$$ | Net Shoreline Movement — total displacement (m) |
| **EPR** | $$\text{NSM} / \Delta t$$ | End Point Rate — displacement per year (m/yr) |
| **LRR** | Least-squares slope across all epochs | Linear Regression Rate (m/yr) with R², p-value, SE |

5. Classifies each transect as accretion / erosion / stable / outlier
6. Exports maps, profiles, and histograms

### Key parameters (editable in script)

| Parameter | Default | Description |
|---|---|---|
| `TRANSECT_SPACING_M` | `500` | Along-shore spacing between transects (m) |
| `TRANSECT_LENGTH_M` | `1000` | Cross-shore length of each transect (m) |
| `OUTLIER_SIGMA` | `2.5` | Z-score threshold for outlier flagging |
| `ACCRETION_MIN_M` | `1.0` | Minimum displacement to classify as accretion (m) |
| `EROSION_MIN_M` | `1.0` | Minimum displacement to classify as erosion (m) |

### Outputs
```
output_dsas/
├── transects.geojson
├── transect_intersections.csv
├── transects_with_metrics.geojson
├── transects_with_metrics.shp
├── transect_origins_metrics.geojson
├── dsas_metrics.csv
├── map_nsm.png
├── map_epr.png
├── map_lrr.png
├── change_profile_all_metrics.png
└── metric_distributions.png
```

---

## Installation

```bash
git clone https://github.com/yourusername/java-coastline-change.git
cd java-coastline-change
pip install -r requirements.txt
```

> **Note:** `pyfes` requires a separate installation step.
> Follow the official guide at https://github.com/CNES/aviso-fes

---

## Running the Pipeline

```bash
# Step 1: Run GEE script in browser, download exports to shorelines/

# Step 2: Tidal correction
cd Part2_TidalCorrection
python tidal_slr_correction.py

# Step 3: DSAS analysis
cd ../Part3_DSAS
python dsas_analysis.py
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `geopandas` | Spatial data I/O and vector operations |
| `shapely` | Geometry construction and intersection |
| `pyfes` | FES2022 tidal constituent evaluation |
| `pyproj` | Geodetic distance and projection |
| `scipy` | Linear regression (LRR) |
| `matplotlib` | All plots and maps |
| `numpy` / `pandas` | Array and tabular data processing |
| `Pillow` | Summary grid image assembly |

---

## References

```
Vos, K., Splinter, K.D., Harley, M.D., Simmons, J.A., Turner, I.L. (2019).
CoastSat: A Google Earth Engine-enabled Python toolkit to extract shorelines
from publicly available satellite imagery.
Environmental Modelling & Software, 122, 104528.
https://doi.org/10.1016/j.envsoft.2019.104528

Himmelstoss, E.A., Henderson, R.E., Kratzmann, M.G., Farris, A.S. (2021).
Digital Shoreline Analysis System (DSAS) Version 5.1 User Guide.
U.S. Geological Survey Open-File Report 2021-1091.
https://doi.org/10.3133/ofr20211091

Lyard, F.H., Allain, D.J., Cancet, M., Carrere, L., Picot, N. (2021).
FES2014 global ocean tides atlas: design and performance.
Ocean Science, 17, 615–649.
https://doi.org/10.5194/os-17-615-2021

Abidin, H.Z., Andreas, H., Gumilar, I., Fukuda, Y., Pohan, Y.E., Deguchi, T. (2011).
Land subsidence of Jakarta (Indonesia) and its relation with urban development.
Natural Hazards, 59(3), 1753–1771.


---

## License

MIT License — free to use and adapt with attribution.
