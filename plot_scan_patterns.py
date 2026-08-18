"""
Plot the 16x16 root-pair and square-pattern scan line CSVs, with the
measurement points overlaid on top.

Requires: pandas, matplotlib
    pip install pandas matplotlib

Notes on the data:
- root_lines.csv / 16x16_lines.csv  -> two columns x,y: the traced scan path.
- root_meas_pts.csv / meas_pts_16x.csv -> two columns, same row count as the
  matching *_lines.csv. Each row is (0,0) except at 256 rows where it is
  (5,5) -- a flag marking "a measurement was taken here". We use those
  flagged row indices to pick out the corresponding (x,y) from *_lines.csv.
"""

import pandas as pd
import matplotlib.pyplot as plt

# --- File paths (edit if you move the files) ---
ROOT_LINES = r"16x16_root_pair\root_lines.csv"
ROOT_MEAS = r"16x16_root_pair\root_meas_pts.csv"

SQUARE_LINES = r"16x16_square_Pattern\16x16_lines.csv"
SQUARE_MEAS = r"16x16_square_Pattern\meas_pts_16x.csv"


def load_pattern(lines_csv, meas_csv):
    lines = pd.read_csv(lines_csv, header=None, names=["x", "y"])
    meas = pd.read_csv(meas_csv, header=None, names=["a", "b"])
    # Flagged rows (a,b) == (5,5) mark measurement points; pull the matching
    # x,y from the line trace at the same row index.
    flagged = meas.index[(meas["a"] != 0) | (meas["b"] != 0)]
    meas_points = lines.loc[flagged]
    return lines, meas_points


root_lines, root_meas = load_pattern(ROOT_LINES, ROOT_MEAS)
square_lines, square_meas = load_pattern(SQUARE_LINES, SQUARE_MEAS)

fig, axes = plt.subplots(1, 2, figsize=(14, 7))

axes[0].plot(root_lines["x"], root_lines["y"], linewidth=0.8, color="tab:blue", zorder=1)
axes[0].scatter(root_meas["x"], root_meas["y"], s=14, color="red", zorder=2, label="Measurement points")
axes[0].set_title(f"Root Pair Pattern ({len(root_meas)} meas. points)")
axes[0].set_xlabel("X")
axes[0].set_ylabel("Y")
axes[0].axis("equal")
axes[0].grid(True, alpha=0.3)
axes[0].legend()

axes[1].plot(square_lines["x"], square_lines["y"], linewidth=0.8, color="tab:orange", zorder=1)
axes[1].scatter(square_meas["x"], square_meas["y"], s=14, color="red", zorder=2, label="Measurement points")
axes[1].set_title(f"16x16 Square Pattern ({len(square_meas)} meas. points)")
axes[1].set_xlabel("X")
axes[1].set_ylabel("Y")
axes[1].axis("equal")
axes[1].grid(True, alpha=0.3)
axes[1].legend()

plt.tight_layout()
plt.savefig("scan_patterns_with_meas_points.png", dpi=150)
plt.show()
