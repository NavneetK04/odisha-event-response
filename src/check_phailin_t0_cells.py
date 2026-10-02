from pathlib import Path
import pickle
import numpy as np
import pandas as pd
import geopandas as gpd
import config as cfg

PROJECT_ROOT = Path(__file__).resolve().parents[1]

impact_file = (
    PROJECT_ROOT
    / "outputs"
    / "impact_benchmark"
    / "phailin_T0_impact.pkl"
)

grid_file = cfg.P1_GRID

# ---------------------------------------------------------
# Load impact
# ---------------------------------------------------------
with open(impact_file, "rb") as f:
    impact = pickle.load(f)

# Impact matrix:
# rows = 500 ensemble members
# columns = 525 exposure/grid cells
loss_matrix = impact.imp_mat

cell_loss_total = np.asarray(
    loss_matrix.sum(axis=0)
).ravel()

cell_loss_mean = cell_loss_total / loss_matrix.shape[0]

# ---------------------------------------------------------
# Load exact P1 grid
# ---------------------------------------------------------
grid = gpd.read_file(grid_file)

grid = grid.sort_values(
    by=["geometry"]
).reset_index(drop=True)

# Use Impact exposure coordinates instead of relying on
# GeoDataFrame row ordering.
coords = np.asarray(impact.coord_exp)

result = pd.DataFrame({
    "lon": coords[:, 1],
    "lat": coords[:, 0],
    "total_loss_all_members": cell_loss_total,
    "mean_loss_per_member": cell_loss_mean,
})

# Keep materially nonzero cells
nonzero = result[
    result["total_loss_all_members"] > 1e3
].copy()

nonzero = nonzero.sort_values(
    "total_loss_all_members",
    ascending=False
)

print("\n========================================")
print("PHAILIN T-0 PER-CELL LOSS CHECK")
print("========================================")

print(f"Exposure cells: {len(result)}")
print(f"Cells with loss > Rs 1,000: {len(nonzero)}")

print("\nNonzero cells:")
print(
    nonzero.to_string(
        index=False,
        formatters={
            "lon": "{:.6f}".format,
            "lat": "{:.6f}".format,
            "total_loss_all_members": "{:,.2f}".format,
            "mean_loss_per_member": "{:,.2f}".format,
        }
    )
)

print("\nTotal across cells:")
print(
    f"Rs {result['total_loss_all_members'].sum():,.2f}"
)

print("\nCell locations:")
for _, r in nonzero.iterrows():
    print(
        f"  ({r.lon:.6f}, {r.lat:.6f}) "
        f"-> mean loss Rs {r.mean_loss_per_member:,.2f}"
    )
