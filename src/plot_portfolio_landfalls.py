from pathlib import Path
import config as cfg
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point


# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

LINKAGE = cfg.P2_LINKAGE

OUT = PROJECT_ROOT / "outputs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------
# 1. Load PRIMARY portfolio cells
# ---------------------------------------------------------
df = pd.read_csv(LINKAGE)

primary = df[df["scenario"].str.upper() == "PRIMARY"].copy()

# Keep only occupied cells
primary = primary[primary["total_tiv"] > 0].copy()

print(f"PRIMARY occupied cells: {len(primary)}")

assert len(primary) == 17, (
    f"Expected 17 occupied PRIMARY cells, found {len(primary)}"
)


# ---------------------------------------------------------
# 2. Landfall points
# ---------------------------------------------------------
landfalls = {
    "Fani": (85.756775, 19.765304),
    "Phailin": (84.954138, 19.288445),
    "Titli": (84.388326, 18.613110),
    "Amphan": (88.198464, 21.814257),
}

lf = pd.DataFrame(
    [
        {"event": event, "lon": lon, "lat": lat}
        for event, (lon, lat) in landfalls.items()
    ]
)

lf_gdf = gpd.GeoDataFrame(
    lf,
    geometry=gpd.points_from_xy(lf.lon, lf.lat),
    crs="EPSG:4326",
)


# ---------------------------------------------------------
# 3. Portfolio cells as points
# ---------------------------------------------------------
portfolio_gdf = gpd.GeoDataFrame(
    primary,
    geometry=gpd.points_from_xy(
        primary["grid_lon"],
        primary["grid_lat"],
    ),
    crs="EPSG:4326",
)


# ---------------------------------------------------------
# 4. Natural Earth coastline
# ---------------------------------------------------------
import cartopy.io.shapereader as shpreader

coast_path = shpreader.natural_earth(
    resolution="10m",
    category="physical",
    name="coastline",
)

coast = gpd.read_file(coast_path)


# ---------------------------------------------------------
# 5. Plot
# ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 9))

coast.plot(
    ax=ax,
    linewidth=0.8,
)

# Portfolio cells
portfolio_gdf.plot(
    ax=ax,
    marker="s",
    markersize=80,
    label="PRIMARY occupied cells (17)",
)

# Landfalls
lf_gdf.plot(
    ax=ax,
    marker="*",
    markersize=180,
    label="Geometric landfall",
)

# Labels
for _, row in lf.iterrows():
    ax.annotate(
        row["event"],
        (row["lon"], row["lat"]),
        xytext=(7, 7),
        textcoords="offset points",
        fontsize=10,
        fontweight="bold",
    )


# ---------------------------------------------------------
# 6. Zoom to Odisha + landfalls
# ---------------------------------------------------------
ax.set_xlim(83.8, 88.6)
ax.set_ylim(17.5, 22.3)

ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")

ax.set_title(
    "Project 2 PRIMARY Portfolio Concentration and Cyclone Landfalls"
)

ax.legend()

ax.grid(True, alpha=0.25)

plt.tight_layout()

outfile = OUT / "portfolio_17_cells_landfalls.png"
plt.savefig(outfile, dpi=300, bbox_inches="tight")

print(f"\nSaved:\n{outfile}")

plt.show()
