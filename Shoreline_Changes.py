import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import unary_union, nearest_points, split, snap
from scipy import stats
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR       = Path(__file__).parent.resolve()
SHORELINE_DIR  = BASE_DIR / "shorelines_corrected"   # output from part 2
OUTPUT_DIR     = BASE_DIR / "output_dsas"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEARS = [1990, 1995, 2002, 2005, 2010, 2015, 2020, 2025]
REF_YEAR       = 1990

TRANSECT_SPACING_M   = 500.0
TRANSECT_LENGTH_M    = 1000.0
TRANSECT_SMOOTH_WIN  = 5
OUTLIER_SIGMA        = 2.5
ACCRETION_MIN_M      = 1.0
EROSION_MIN_M        = 1.0

EPOCH_COLORS = {
    1990: "#1abc9c", 1995: "#16a085", 2002: "#2c3e50",
    2005: "#2980b9", 2010: "#27ae60", 2015: "#f39c12",
    2020: "#e67e22", 2025: "#8e44ad",
}
CHANGE_COLORS = {
    "accretion": "#2ecc71",
    "erosion":   "#e74c3c",
    "stable":    "#95a5a6",
    "outlier":   "#f39c12",
}

# =============================================================================
# HELPERS
# =============================================================================

def meters_to_degrees(meters, lat=0.0):
    return meters / (111320 * np.cos(np.radians(lat)))

def classify_change(delta_m, z_score):
    if abs(z_score) > OUTLIER_SIGMA:
        return "outlier"
    elif delta_m >= ACCRETION_MIN_M:
        return "accretion"
    elif delta_m <= -EROSION_MIN_M:
        return "erosion"
    return "stable"

def smooth_line(coords, window=5):
    if len(coords) < window:
        return coords
    coords = np.array(coords)
    half   = window // 2
    smoothed = []
    for i in range(len(coords)):
        lo = max(0, i - half)
        hi = min(len(coords), i + half + 1)
        smoothed.append(coords[lo:hi].mean(axis=0))
    return smoothed

# =============================================================================
# LOAD SHORELINES
# =============================================================================

def load_shorelines():
    shorelines = {}
    print(f"\n{'='*65}")
    print(f"  LOADING SHORELINES")
    print(f"{'='*65}")
    for year in YEARS:
        candidates = (list(SHORELINE_DIR.glob(f"*{year}*.shp")) +
                      list(SHORELINE_DIR.glob(f"*{year}*.geojson")))
        if not candidates:
            print(f"   WARNING: No file found for {year}, skipping.")
            continue
        gdf = gpd.read_file(candidates[0])
        if gdf.crs is None:
            gdf = gdf.set_crs("EPSG:4326")
        elif gdf.crs.to_epsg() != 4326:
            gdf = gdf.to_crs("EPSG:4326")
        gdf = gdf[gdf.geometry.notnull()].copy()
        shorelines[year] = gdf
        print(f"   {year}: {candidates[0].name}  ({len(gdf)} features)")
    return shorelines

# =============================================================================
# BUILD BASELINE FROM REFERENCE SHORELINE
# =============================================================================

def build_baseline(shorelines):
    print(f"\n  Building baseline from {REF_YEAR} shoreline...")
    ref = shorelines[REF_YEAR]
    geoms = ref.geometry.tolist()
    merged = unary_union(geoms)
    if merged.geom_type == "MultiLineString":
        lines = list(merged.geoms)
    elif merged.geom_type == "LineString":
        lines = [merged]
    else:
        lines = [g for g in merged.geoms if g.geom_type == "LineString"]
    print(f"   Baseline: {len(lines)} segment(s)")
    return lines

# =============================================================================
# GENERATE TRANSECTS
# =============================================================================

def generate_transects(baseline_lines):
    print(f"\n  Generating transects...")
    print(f"   Spacing: {TRANSECT_SPACING_M} m  |  Length: {TRANSECT_LENGTH_M} m")

    transects  = []
    transect_id = 0

    for line in baseline_lines:
        coords = list(line.coords)
        if len(coords) < 2:
            continue

        smoothed = smooth_line(coords, TRANSECT_SMOOTH_WIN)
        total_len = line.length

        dist = 0.0
        while dist <= total_len:
            pt = line.interpolate(dist)
            lat = pt.y

            spacing_deg = meters_to_degrees(TRANSECT_SPACING_M, lat)
            length_deg  = meters_to_degrees(TRANSECT_LENGTH_M,  lat)

            frac_lo = max(0.0, dist - spacing_deg * 0.01)
            frac_hi = min(total_len, dist + spacing_deg * 0.01)
            pt_lo   = line.interpolate(frac_lo)
            pt_hi   = line.interpolate(frac_hi)

            dx = pt_hi.x - pt_lo.x
            dy = pt_hi.y - pt_lo.y
            norm = np.sqrt(dx**2 + dy**2)

            if norm < 1e-10:
                dist += meters_to_degrees(TRANSECT_SPACING_M, lat)
                continue

            perp_x = -dy / norm
            perp_y =  dx / norm

            half = length_deg / 2.0
            p1   = Point(pt.x + perp_x * half, pt.y + perp_y * half)
            p2   = Point(pt.x - perp_x * half, pt.y - perp_y * half)

            transects.append({
                "transect_id": transect_id,
                "origin_x":    pt.x,
                "origin_y":    pt.y,
                "geometry":    LineString([p1, p2]),
            })
            transect_id += 1
            dist += meters_to_degrees(TRANSECT_SPACING_M, lat)

    gdf = gpd.GeoDataFrame(transects, crs="EPSG:4326")
    print(f"   Generated {len(gdf)} transects")
    return gdf

# =============================================================================
# INTERSECT TRANSECTS WITH SHORELINES
# =============================================================================

def intersect_transects(transects_gdf, shorelines):
    print(f"\n  Intersecting transects with shorelines...")

    results = {tid: {"transect_id": tid,
                     "origin_x": row["origin_x"],
                     "origin_y": row["origin_y"]}
               for tid, row in transects_gdf.iterrows()}

    for year, sl_gdf in shorelines.items():
        sl_union = unary_union(sl_gdf.geometry.tolist())
        hit = 0

        for tid, row in transects_gdf.iterrows():
            transect = row["geometry"]
            origin   = Point(row["origin_x"], row["origin_y"])

            try:
                intersection = transect.intersection(sl_union)
            except Exception:
                results[tid][f"dist_{year}_m"] = np.nan
                continue

            if intersection.is_empty:
                results[tid][f"dist_{year}_m"] = np.nan
                continue

            if intersection.geom_type == "Point":
                pts = [intersection]
            elif intersection.geom_type == "MultiPoint":
                pts = list(intersection.geoms)
            elif intersection.geom_type in ("LineString", "MultiLineString",
                                             "GeometryCollection"):
                pts = []
                geoms = (list(intersection.geoms)
                         if hasattr(intersection, "geoms") else [intersection])
                for g in geoms:
                    if g.geom_type == "Point":
                        pts.append(g)
                    elif hasattr(g, "coords"):
                        pts += [Point(c) for c in g.coords]
            else:
                results[tid][f"dist_{year}_m"] = np.nan
                continue

            if not pts:
                results[tid][f"dist_{year}_m"] = np.nan
                continue

            nearest = min(pts, key=lambda p: origin.distance(p))
            lat_deg = origin.y
            dist_m  = origin.distance(nearest) * 111320 * np.cos(np.radians(lat_deg))

            side_vec = np.array([nearest.x - origin.x, nearest.y - origin.y])
            tcoords  = np.array(transect.coords)
            t_vec    = tcoords[-1] - tcoords[0]
            sign     = 1.0 if np.dot(side_vec, t_vec) >= 0 else -1.0

            results[tid][f"dist_{year}_m"]    = float(sign * dist_m)
            results[tid][f"int_x_{year}"]     = float(nearest.x)
            results[tid][f"int_y_{year}"]     = float(nearest.y)
            hit += 1

        print(f"   {year}: {hit}/{len(transects_gdf)} transects intersected")

    return pd.DataFrame(list(results.values()))

# =============================================================================
# COMPUTE DSAS METRICS: EPR, LRR, NSM
# =============================================================================

def compute_dsas_metrics(df):
    print(f"\n  Computing DSAS metrics (EPR, LRR, NSM)...")

    available_years = [y for y in YEARS if f"dist_{y}_m" in df.columns]
    print(f"   Years with data: {available_years}")

    records = []

    for _, row in df.iterrows():
        tid    = row["transect_id"]
        ox, oy = row["origin_x"], row["origin_y"]

        year_vals = {}
        for y in available_years:
            v = row.get(f"dist_{y}_m", np.nan)
            if not np.isnan(v):
                year_vals[y] = v

        if len(year_vals) < 2:
            continue

        yrs  = np.array(sorted(year_vals.keys()), dtype=float)
        dsts = np.array([year_vals[y] for y in yrs.astype(int)])

        # NSM — Net Shoreline Movement (first to last available year)
        nsm_m   = dsts[-1] - dsts[0]
        nsm_yrs = yrs[-1] - yrs[0]

        # EPR — End Point Rate (m/yr)
        epr = nsm_m / nsm_yrs if nsm_yrs > 0 else np.nan

        # LRR — Linear Regression Rate (m/yr)
        if len(yrs) >= 3:
            slope, intercept, r_val, p_val, std_err = stats.linregress(yrs, dsts)
            lrr     = float(slope)
            lrr_r2  = float(r_val**2)
            lrr_p   = float(p_val)
            lrr_se  = float(std_err)
        else:
            lrr, lrr_r2, lrr_p, lrr_se = np.nan, np.nan, np.nan, np.nan

        rec = {
            "transect_id": tid,
            "origin_x":    ox,
            "origin_y":    oy,
            "nsm_m":       float(nsm_m),
            "nsm_years":   float(nsm_yrs),
            "epr_m_yr":    float(epr) if not np.isnan(epr) else np.nan,
            "lrr_m_yr":    lrr,
            "lrr_r2":      lrr_r2,
            "lrr_p":       lrr_p,
            "lrr_se":      lrr_se,
            "n_years":     len(year_vals),
        }

        for y in available_years:
            rec[f"dist_{y}_m"] = year_vals.get(y, np.nan)

        records.append(rec)

    df_metrics = pd.DataFrame(records)

    for col in ["nsm_m", "epr_m_yr", "lrr_m_yr"]:
        std = df_metrics[col].std()
        if std > 0:
            df_metrics[f"z_{col}"] = (df_metrics[col] - df_metrics[col].mean()) / std
        else:
            df_metrics[f"z_{col}"] = 0.0

    df_metrics["nsm_class"] = df_metrics.apply(
        lambda r: classify_change(r["nsm_m"], r["z_nsm_m"]), axis=1)
    df_metrics["epr_class"] = df_metrics.apply(
        lambda r: classify_change(r["epr_m_yr"], r.get("z_epr_m_yr", 0)), axis=1)
    df_metrics["lrr_class"] = df_metrics.apply(
        lambda r: classify_change(r["lrr_m_yr"], r.get("z_lrr_m_yr", 0))
        if not np.isnan(r["lrr_m_yr"]) else "stable", axis=1)

    print(f"   Metrics computed for {len(df_metrics)} transects")
    print(f"\n── Summary ──────────────────────────────────────────────────")
    for metric, col, cls_col in [
        ("NSM",  "nsm_m",    "nsm_class"),
        ("EPR",  "epr_m_yr", "epr_class"),
        ("LRR",  "lrr_m_yr", "lrr_class"),
    ]:
        v = df_metrics[col].dropna()
        c = df_metrics[cls_col].value_counts()
        print(f"   {metric}: mean={v.mean():+.2f}  std={v.std():.2f}  "
              f"accretion={c.get('accretion',0)}  "
              f"erosion={c.get('erosion',0)}  "
              f"stable={c.get('stable',0)}  "
              f"outlier={c.get('outlier',0)}")
    print(f"─────────────────────────────────────────────────────────────")

    return df_metrics

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

def save_outputs(df_metrics, transects_gdf):
    df_metrics.to_csv(OUTPUT_DIR / "dsas_metrics.csv", index=False)
    print(f"\n   Saved dsas_metrics.csv  ({len(df_metrics)} transects)")

    merged = transects_gdf.merge(df_metrics, on="transect_id", how="left")
    merged.to_file(OUTPUT_DIR / "transects_with_metrics.geojson", driver="GeoJSON")
    merged.to_file(OUTPUT_DIR / "transects_with_metrics.shp")
    print(f"   Saved transects_with_metrics.geojson / .shp")

    origins = gpd.GeoDataFrame(
        df_metrics,
        geometry=[Point(r["origin_x"], r["origin_y"]) for _, r in df_metrics.iterrows()],
        crs="EPSG:4326"
    )
    origins.to_file(OUTPUT_DIR / "transect_origins_metrics.geojson", driver="GeoJSON")
    print(f"   Saved transect_origins_metrics.geojson")

# =============================================================================
# PLOTS
# =============================================================================

def plot_metric_map(df_metrics, metric_col, class_col, title, filename):
    fig, ax = plt.subplots(figsize=(16, 10))
    fig.patch.set_facecolor('#FAFAFA')
    ax.set_facecolor('#EEF2F7')

    valid = df_metrics[df_metrics[metric_col].notna()].copy()
    vmin  = valid[metric_col].quantile(0.05)
    vmax  = valid[metric_col].quantile(0.95)
    norm  = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    sc = ax.scatter(valid["origin_x"], valid["origin_y"],
                    c=valid[metric_col], cmap="RdYlGn",
                    norm=norm, s=18, alpha=0.85,
                    edgecolors='none', zorder=4)

    cbar = plt.colorbar(sc, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label(f"{metric_col} [m or m/yr]", fontsize=10)

    for cls, color in CHANGE_COLORS.items():
        mask = valid[class_col] == cls
        if mask.sum() == 0:
            continue
        ax.scatter([], [], c=color, s=40, label=f"{cls.capitalize()} ({mask.sum()})")

    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
    ax.set_xlabel("Longitude (°)", fontsize=11)
    ax.set_ylabel("Latitude (°)", fontsize=11)
    ax.legend(fontsize=9, framealpha=0.9, loc='upper right')
    ax.grid(True, alpha=0.2, linestyle='--')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved {filename}")


def plot_change_profile(df_metrics):
    available_years = [y for y in YEARS if f"dist_{y}_m" in df_metrics.columns]
    x = np.arange(len(df_metrics))

    fig, axes = plt.subplots(3, 1, figsize=(18, 14))
    fig.patch.set_facecolor('#FAFAFA')

    ax1 = axes[0]
    ax1.set_facecolor('#F8F8F8')
    for year in available_years:
        col = f"dist_{year}_m"
        vals = df_metrics[col].values
        ax1.plot(x, vals, color=EPOCH_COLORS.get(year, "#333333"),
                 linewidth=1.2, alpha=0.8, label=str(year))
    ax1.axhline(0, color='black', linewidth=0.8)
    ax1.set_title("Shoreline Position per Epoch (Distance from Baseline)", fontsize=12, fontweight='bold')
    ax1.set_ylabel("Distance [m]", fontsize=10)
    ax1.legend(fontsize=8, ncol=4, framealpha=0.9)
    ax1.grid(True, alpha=0.2, linestyle='--')

    ax2 = axes[1]
    ax2.set_facecolor('#F8F8F8')
    nsm  = df_metrics["nsm_m"].values
    cols = [CHANGE_COLORS.get(c, "#aaaaaa") for c in df_metrics["nsm_class"]]
    ax2.bar(x, nsm, color=cols, alpha=0.80, edgecolor='none', width=0.85)
    ax2.axhline(0, color='black', linewidth=0.8)
    ax2.axhline(nsm[~np.isnan(nsm)].mean(), color='navy', linewidth=1.5,
                linestyle='--', label=f"Mean NSM: {np.nanmean(nsm):+.1f} m")
    ax2.set_title(f"Net Shoreline Movement (NSM)  |  {available_years[0]}→{available_years[-1]}", fontsize=12, fontweight='bold')
    ax2.set_ylabel("NSM [m]", fontsize=10)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.2, linestyle='--', axis='y')
    patches = [mpatches.Patch(color=v, label=k.capitalize()) for k, v in CHANGE_COLORS.items()]
    ax2.legend(handles=patches, fontsize=8, framealpha=0.9)

    ax3 = axes[2]
    ax3.set_facecolor('#F8F8F8')
    epr  = df_metrics["epr_m_yr"].values
    lrr  = df_metrics["lrr_m_yr"].values
    ax3.plot(x, epr, color="#2980b9", linewidth=1.2, alpha=0.85, label=f"EPR (mean={np.nanmean(epr):+.2f} m/yr)")
    ax3.plot(x, lrr, color="#8e44ad", linewidth=1.2, alpha=0.85, label=f"LRR (mean={np.nanmean(lrr):+.2f} m/yr)")
    ax3.axhline(0, color='black', linewidth=0.8)
    ax3.fill_between(x, epr, 0, where=epr >= 0, alpha=0.10, color="#2ecc71")
    ax3.fill_between(x, epr, 0, where=epr <  0, alpha=0.10, color="#e74c3c")
    ax3.set_title("EPR & LRR per Transect", fontsize=12, fontweight='bold')
    ax3.set_xlabel("Transect Index", fontsize=10)
    ax3.set_ylabel("Rate [m/yr]", fontsize=10)
    ax3.legend(fontsize=9, framealpha=0.9)
    ax3.grid(True, alpha=0.2, linestyle='--')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "change_profile_all_metrics.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved change_profile_all_metrics.png")


def plot_lrr_histogram(df_metrics):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.patch.set_facecolor('#FAFAFA')

    for ax, (col, label, cls_col) in zip(axes, [
        ("nsm_m",    "NSM [m]",     "nsm_class"),
        ("epr_m_yr", "EPR [m/yr]",  "epr_class"),
        ("lrr_m_yr", "LRR [m/yr]",  "lrr_class"),
    ]):
        vals = df_metrics[col].dropna()
        ax.set_facecolor('#F8F8F8')
        ax.hist(vals[vals >= 0], bins=30, color=CHANGE_COLORS["accretion"],
                alpha=0.75, label="Accretion", edgecolor='white')
        ax.hist(vals[vals <  0], bins=30, color=CHANGE_COLORS["erosion"],
                alpha=0.75, label="Erosion",   edgecolor='white')
        ax.axvline(vals.mean(),   color='navy',   linewidth=1.5, linestyle='--',
                   label=f"Mean: {vals.mean():+.2f}")
        ax.axvline(vals.median(), color='orange', linewidth=1.5, linestyle=':',
                   label=f"Median: {vals.median():+.2f}")
        ax.set_title(f"{label} Distribution", fontsize=11, fontweight='bold')
        ax.set_xlabel(label, fontsize=10)
        ax.set_ylabel("Count", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2, linestyle='--', axis='y')

    plt.suptitle("Shoreline Change Metric Distributions  |  Java Coast 1990–2025",
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "metric_distributions.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved metric_distributions.png")

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 65)
    print("  Java Coastline Change 1990-2025")
    print("  DSAS-Equivalent: EPR / LRR / NSM Analysis")
    print("=" * 65)

    shorelines   = load_shorelines()
    if not shorelines:
        raise RuntimeError("No shoreline files found. Check SHORELINE_DIR.")

    baseline_lines = build_baseline(shorelines)
    transects_gdf  = generate_transects(baseline_lines)
    transects_gdf.to_file(OUTPUT_DIR / "transects.geojson", driver="GeoJSON")
    print(f"   Saved transects.geojson")

    df_intersect = intersect_transects(transects_gdf, shorelines)
    df_intersect.to_csv(OUTPUT_DIR / "transect_intersections.csv", index=False)
    print(f"   Saved transect_intersections.csv")

    df_metrics = compute_dsas_metrics(df_intersect)
    save_outputs(df_metrics, transects_gdf)

    print(f"\n  Generating plots...")
    plot_metric_map(df_metrics, "nsm_m",    "nsm_class",
                    "Net Shoreline Movement (NSM)  |  Java 1990–2025",
                    "map_nsm.png")
    plot_metric_map(df_metrics, "epr_m_yr", "epr_class",
                    "End Point Rate (EPR)  |  Java 1990–2025",
                    "map_epr.png")
    plot_metric_map(df_metrics, "lrr_m_yr", "lrr_class",
                    "Linear Regression Rate (LRR)  |  Java 1990–2025",
                    "map_lrr.png")
    plot_change_profile(df_metrics)
    plot_lrr_histogram(df_metrics)

    print(f"\n{'='*65}")
    print(f"  DONE -> {OUTPUT_DIR}")
    print(f"{'='*65}")
    print(f"\n  Outputs:")
    print(f"    transects.geojson")
    print(f"    transect_intersections.csv")
    print(f"    transects_with_metrics.geojson / .shp")
    print(f"    transect_origins_metrics.geojson")
    print(f"    dsas_metrics.csv")
    print(f"    map_nsm.png / map_epr.png / map_lrr.png")
    print(f"    change_profile_all_metrics.png")
    print(f"    metric_distributions.png")

if __name__ == "__main__":
    main()
