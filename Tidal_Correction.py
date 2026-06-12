import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from shapely.geometry import Point, LineString
from datetime import datetime, timedelta
from pathlib import Path
from pyproj import Geod
import warnings

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR        = Path(__file__).parent.resolve()
SHORELINE_DIR   = BASE_DIR / "shorelines"
FES_YAML        = BASE_DIR / "fes2022_indonesia_clipped.yaml"
OUTPUT_DIR      = BASE_DIR / "output_tidal_corrected"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TIDAL_WINDOWS = {
    1990: ("1989-01-01", "1990-12-31"),
    1995: ("1994-01-01", "1995-12-31"),
    2002: ("2001-01-01", "2002-12-31"),
    2005: ("2004-01-01", "2005-12-31"),
    2010: ("2009-01-01", "2010-12-31"),
    2015: ("2014-01-01", "2015-12-31"),
    2020: ("2019-01-01", "2020-12-31"),
    2025: ("2024-01-01", "2025-12-31"),
}

SENSOR_CONFIG = {
    1990: {"sensor": "Landsat-4/5", "pixel_m": 30.0},
    1995: {"sensor": "Landsat-4/5", "pixel_m": 30.0},
    2002: {"sensor": "Landsat-5/7", "pixel_m": 30.0},
    2005: {"sensor": "Landsat-5/7", "pixel_m": 30.0},
    2010: {"sensor": "Landsat-5/7", "pixel_m": 30.0},
    2015: {"sensor": "Landsat-8",   "pixel_m": 30.0},
    2020: {"sensor": "Landsat-8",   "pixel_m": 30.0},
    2025: {"sensor": "Landsat-8/9", "pixel_m": 30.0},
}

SLR_RATE_MM_YR  = 3.7
SUBSIDENCE_RATE_MM_YR = {
    "micro-tidal": 4.5,
    "meso-tidal":  12.0,
    "macro-tidal": 7.5,
}
REF_YEAR = 2002

EPOCH_COLORS = {
    1990: "#1abc9c", 1995: "#16a085", 2002: "#2c3e50",
    2005: "#2980b9", 2010: "#27ae60", 2015: "#f39c12",
    2020: "#e67e22", 2025: "#8e44ad",
}

# =============================================================================
# HELPERS
# =============================================================================

_FES_CONFIG_CACHE = None

def resolve_fes_yaml(yaml_path):
    yaml_dir = yaml_path.parent.resolve()
    resolved = yaml_dir / "fes2022_resolved.yaml"
    content  = yaml_path.read_text()
    resolved.write_text(content.replace("./", str(yaml_dir).replace("\\", "/") + "/"))
    return resolved

def get_fes_config(resolved_yaml):
    global _FES_CONFIG_CACHE
    if _FES_CONFIG_CACHE is None:
        import pyfes
        _FES_CONFIG_CACHE = pyfes.config.load(str(resolved_yaml))
    return _FES_CONFIG_CACHE

def tidal_class(tidal_range_m):
    if tidal_range_m < 2.0:
        return "micro-tidal"
    elif tidal_range_m < 4.0:
        return "meso-tidal"
    return "macro-tidal"

def compute_tidal_stats(lon, lat, date_start, date_end, resolved_yaml, config=None):
    import pyfes
    t_start = datetime.strptime(date_start, "%Y-%m-%d")
    t_end   = datetime.strptime(date_end,   "%Y-%m-%d")
    timestamps = []
    t = t_start
    while t <= t_end:
        timestamps.append(t)
        t += timedelta(hours=1)
    n        = len(timestamps)
    times_np = np.array([np.datetime64(ts.strftime("%Y-%m-%dT%H:%M:%S")) for ts in timestamps], dtype="datetime64[s]")
    lons     = np.full(n, lon, dtype=np.float64)
    lats     = np.full(n, lat, dtype=np.float64)
    if config is None:
        config = get_fes_config(resolved_yaml)
    ocean_sp, lp, _ = pyfes.evaluate_tide(config["tide"],   times_np, lons, lats)
    load_sp,  _,  _ = pyfes.evaluate_tide(config["radial"], times_np, lons, lats)
    total_m = (ocean_sp + load_sp) / 100.0
    valid   = total_m[~np.isnan(total_m)]
    if len(valid) == 0:
        return {"hat_m": 0.0, "lat_m": 0.0, "msl_m": 0.0, "range_m": 0.0, "std_m": 0.0, "n_samples": 0}
    return {
        "hat_m":     float(np.max(valid)),
        "lat_m":     float(np.min(valid)),
        "msl_m":     float(np.mean(valid)),
        "range_m":   float(np.max(valid) - np.min(valid)),
        "std_m":     float(np.std(valid)),
        "n_samples": int(len(valid)),
    }

def compute_vertical_offset(msl_epoch, msl_ref, year, ref_year=REF_YEAR):
    dt_years     = year - ref_year
    slr_m        = (SLR_RATE_MM_YR * dt_years) / 1000.0
    tidal_delta  = msl_epoch - msl_ref
    return float(tidal_delta + slr_m)

def subsidence_m(tidal_cls, year, ref_year=REF_YEAR):
    dt   = year - ref_year
    rate = SUBSIDENCE_RATE_MM_YR.get(tidal_cls, 8.0)
    return (rate * dt) / 1000.0

def lwl_from_stats(stats):
    return float(stats["lat_m"])

# =============================================================================
# LOAD GEE SHORELINES
# =============================================================================

def load_shorelines():
    shorelines = {}
    for year in TIDAL_WINDOWS:
        candidates = list(SHORELINE_DIR.glob(f"*{year}*.shp")) + \
                     list(SHORELINE_DIR.glob(f"*{year}*.geojson"))
        if not candidates:
            print(f"   WARNING: No shoreline file found for {year}, skipping.")
            continue
        gdf = gpd.read_file(candidates[0])
        if gdf.crs is None:
            gdf = gdf.set_crs("EPSG:4326")
        elif gdf.crs.to_epsg() != 4326:
            gdf = gdf.to_crs("EPSG:4326")
        shorelines[year] = gdf
        print(f"   Loaded {year}: {candidates[0].name}  ({len(gdf)} features)")
    return shorelines

# =============================================================================
# SAMPLE TIDAL STATS ALONG SHORELINE
# =============================================================================

def sample_tidal_stats_along_shoreline(shorelines, resolved_yaml):
    fes_config = get_fes_config(resolved_yaml)
    results    = {}

    for year, gdf in shorelines.items():
        d_start, d_end = TIDAL_WINDOWS[year]
        print(f"\n  [{year}] Sampling tidal stats along {len(gdf)} features...")
        records = []

        for idx, row in gdf.iterrows():
            geom = row.geometry
            if geom is None:
                continue
            centroid = geom.centroid
            lon, lat = centroid.x, centroid.y

            try:
                stats = compute_tidal_stats(lon, lat, d_start, d_end, resolved_yaml, config=fes_config)
            except Exception as e:
                print(f"     WARNING idx={idx}: {e}")
                stats = {"hat_m": 0.0, "lat_m": 0.0, "msl_m": 0.0, "range_m": 0.0, "std_m": 0.0, "n_samples": 0}

            tc    = tidal_class(stats["range_m"])
            subs  = subsidence_m(tc, year)
            lwl   = lwl_from_stats(stats)

            records.append({
                "idx":          idx,
                "lon":          lon,
                "lat":          lat,
                "year":         year,
                "sensor":       SENSOR_CONFIG.get(year, {}).get("sensor", "Unknown"),
                "pixel_m":      SENSOR_CONFIG.get(year, {}).get("pixel_m", 30.0),
                "hat_m":        stats["hat_m"],
                "lat_m":        stats["lat_m"],
                "msl_m":        stats["msl_m"],
                "range_m":      stats["range_m"],
                "std_m":        stats["std_m"],
                "tidal_class":  tc,
                "lwl_m":        lwl,
                "subsidence_m": subs,
                "n_samples":    stats["n_samples"],
                "geometry":     geom,
            })

        results[year] = pd.DataFrame(records)
        print(f"   Done {year}: {len(records)} records")

    return results

# =============================================================================
# APPLY TIDAL + SLR VERTICAL OFFSET CORRECTION
# =============================================================================

def apply_tidal_slr_correction(tidal_results):
    ref_msl_lookup = {}

    if REF_YEAR in tidal_results:
        ref_df = tidal_results[REF_YEAR]
        ref_msl_lookup = dict(zip(ref_df["idx"], ref_df["msl_m"]))
        ref_msl_global = float(ref_df["msl_m"].mean())
    else:
        ref_msl_global = 0.0
        print(f"   WARNING: REF_YEAR {REF_YEAR} not in tidal_results. Using global mean=0.")

    corrected = {}

    for year, df in tidal_results.items():
        df = df.copy()
        vertical_offsets = []
        lwl_corrected    = []

        for _, row in df.iterrows():
            ref_msl = ref_msl_lookup.get(row["idx"], ref_msl_global)
            v_offset = compute_vertical_offset(row["msl_m"], ref_msl, year)
            lwl_corr = row["lwl_m"] - v_offset - row["subsidence_m"]
            vertical_offsets.append(v_offset)
            lwl_corrected.append(lwl_corr)

        df["vertical_offset_m"] = vertical_offsets
        df["lwl_corrected_m"]   = lwl_corrected
        df["correction_note"]   = (
            f"LWL = FES2022_LAT - vertical_offset(tidal+SLR@{SLR_RATE_MM_YR}mm/yr) - subsidence"
        )
        corrected[year] = df

    return corrected

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

def save_outputs(corrected, shorelines):
    all_records = []

    for year, df in corrected.items():
        gdf = gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")
        gdf.to_file(OUTPUT_DIR / f"shoreline_tidal_corrected_{year}.geojson", driver="GeoJSON")
        gdf.to_file(OUTPUT_DIR / f"shoreline_tidal_corrected_{year}.shp")
        all_records.append(df.drop(columns=["geometry"]))
        print(f"   Saved {year}: shoreline_tidal_corrected_{year}.geojson / .shp")

    pd.concat(all_records, ignore_index=True).to_csv(
        OUTPUT_DIR / "tidal_correction_all_years.csv", index=False)
    print(f"   Saved tidal_correction_all_years.csv")

# =============================================================================
# PLOT OVERVIEW
# =============================================================================

def plot_overview(corrected):
    years = sorted(corrected.keys())

    fig, axes = plt.subplots(len(years), 1, figsize=(16, 4 * len(years)), sharex=False)
    if len(years) == 1:
        axes = [axes]
    fig.patch.set_facecolor('#FAFAFA')
    fig.suptitle("FES2022 Tidal LWL Correction per Year\n(Tidal ΔMSL + SLR + Subsidence)", fontsize=13, fontweight='bold')

    for ax, year in zip(axes, years):
        df   = corrected[year].reset_index(drop=True)
        x    = np.arange(len(df))
        vals = df["lwl_corrected_m"].values
        ax.bar(x, vals, color=EPOCH_COLORS.get(year, "#333333"), alpha=0.75, edgecolor='none', width=0.85)
        ax.axhline(vals.mean(), color='black', linewidth=1.5, linestyle='--',
                   label=f"Mean LWL corrected: {vals.mean():+.4f} m")
        ax.set_title(f"{year}  |  {SENSOR_CONFIG.get(year,{}).get('sensor','')}  |  n={len(df)}", fontsize=11, fontweight='bold')
        ax.set_ylabel("LWL Corrected [m]", fontsize=10)
        ax.set_xlabel("Feature Index", fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.25, linestyle='--', axis='y')
        ax.set_facecolor('#F8F8F8')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "lwl_correction_overview.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved lwl_correction_overview.png")

    fig, ax = plt.subplots(figsize=(16, 6))
    fig.patch.set_facecolor('#FAFAFA')
    ax.set_facecolor('#F8F8F8')
    for year in years:
        df = corrected[year].reset_index(drop=True)
        ax.plot(np.arange(len(df)), df["lwl_corrected_m"],
                color=EPOCH_COLORS.get(year, "#333333"), linewidth=1.5, alpha=0.85,
                label=f"{year} (mean={df['lwl_corrected_m'].mean():+.4f}m)")
    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_title("LWL Corrected — All Years Comparison\nFES2022 + SLR + Subsidence", fontsize=13, fontweight='bold')
    ax.set_xlabel("Feature Index", fontsize=11)
    ax.set_ylabel("LWL Corrected [m]", fontsize=11)
    ax.legend(fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.25, linestyle='--')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "lwl_all_years_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved lwl_all_years_comparison.png")

# =============================================================================
# PRINT SUMMARY
# =============================================================================

def print_summary(corrected):
    print(f"\n{'='*70}")
    print(f"  TIDAL + SLR CORRECTION SUMMARY")
    print(f"{'='*70}")
    for year in sorted(corrected.keys()):
        df = corrected[year]
        print(f"\n  {year}  [{SENSOR_CONFIG.get(year,{}).get('sensor','')}]:")
        print(f"     Features         : {len(df)}")
        print(f"     Mean MSL         : {df['msl_m'].mean():+.5f} m")
        print(f"     Mean LWL (LAT)   : {df['lat_m'].mean():+.5f} m")
        print(f"     Mean V.Offset    : {df['vertical_offset_m'].mean():+.5f} m")
        print(f"     Mean Subsidence  : {df['subsidence_m'].mean():+.5f} m")
        print(f"     Mean LWL Corr.   : {df['lwl_corrected_m'].mean():+.5f} m")
        print(f"     Tidal classes    : {df['tidal_class'].value_counts().to_dict()}")
    print(f"{'='*70}")

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("  Java Coastline Change 1990-2025")
    print("  Tidal + SLR Vertical Offset Correction  |  FES2022")
    print("=" * 70)

    resolved_yaml = resolve_fes_yaml(FES_YAML)
    shorelines    = load_shorelines()

    if not shorelines:
        raise RuntimeError("No shoreline files found. Check SHORELINE_DIR.")

    tidal_results = sample_tidal_stats_along_shoreline(shorelines, resolved_yaml)
    corrected     = apply_tidal_slr_correction(tidal_results)

    save_outputs(corrected, shorelines)
    print_summary(corrected)
    plot_overview(corrected)

    print(f"\n  DONE -> {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
