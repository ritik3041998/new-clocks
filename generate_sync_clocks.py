"""
generate_sync_clocks.py

Generates pixel / line / frame synchronization clocks for an NLOS galvo
scan, from the same sample-indexed tables that already describe the scan
trajectory (lines.csv) and the measurement/pixel locations (meas_pts.csv).

Core idea (see write-up): pixel/line/frame clocks are NOT derived from an
assumed constant velocity or a fixed sample-per-pixel counter. They are
looked up, sample-by-sample, from the same index space that the geometric
scan-pattern correction already produced. This works identically for a
uniform square-voltage scan and a nonuniform corrected-trapezoidal scan --
only the table content differs, not the logic.

Input file format (matches your existing data): two columns, no header.
    lines.csv     -> x, y                  (galvo voltage sample per row)
    meas_pts.csv  -> a, b                  (0,0 = not a pixel; nonzero = pixel)
Both files must have the same number of rows as the lines.csv they pair with
(row n in meas_pts.csv corresponds to row n / sample n in lines.csv).

Usage
-----
    python generate_sync_clocks.py \
        --lines 16x16_root_pair/root_lines.csv \
        --meas 16x16_root_pair/root_meas_pts.csv \
        --fs 1000 \
        --out out_root

    python generate_sync_clocks.py \
        --lines 16x16_square_Pattern/16x16_lines.csv \
        --meas 16x16_square_Pattern/meas_pts_16x.csv \
        --fs 1000 \
        --out out_square

Run with no arguments to process both of your existing 16x16 datasets
(root pair + square pattern) at 1000 Hz, using their known folder layout.

Outputs (written into --out directory)
---------------------------------------
    sync_clocks_full.csv   one row per DAC sample: n, t_s, x, y,
                            pixel_clk, pixel_index, line_clk, line_index,
                            frame_clk, pixel_valid_gate
    sync_pixel_edges.csv   one row per pixel only: pixel_index, sample_n,
                            t_s, x, y, line_index, is_line_start
    sync_clocks.png        plot of pixel/line/frame clock vs. time, for a
                            visual sanity check
    summary.txt            frame duration, pixel/line counts, dwell-time
                            and gap statistics, at the chosen fs
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_pair(lines_csv, meas_csv):
    lines = pd.read_csv(lines_csv, header=None, names=["x", "y"])
    meas = pd.read_csv(meas_csv, header=None, names=["a", "b"])
    if len(lines) != len(meas):
        raise ValueError(
            f"Row count mismatch: {lines_csv} has {len(lines)} rows, "
            f"{meas_csv} has {len(meas)} rows. They must be sample-aligned."
        )
    return lines, meas


def detect_line_boundaries(flagged_idx, points_per_line=None, gap_factor=1.3):
    """
    Decide which pixels start a new scan line.

    If points_per_line is given, lines are simply chunks of that length
    (matches a known N x N raster where every line has the same pixel count).

    Otherwise, line boundaries are inferred from the sample-index GAP between
    consecutive pixels: within a line, consecutive pixels are close together
    in sample index; at a line turnaround (flyback/settle), the gap is
    noticeably larger. A gap larger than gap_factor * median(gap) marks the
    start of a new line. This generalizes to nonuniform (trapezoidal)
    spacing, since it only compares each gap to the typical *intra-line*
    gap, not to an assumed constant.
    """
    n = len(flagged_idx)
    line_start = np.zeros(n, dtype=bool)
    line_start[0] = True

    if points_per_line is not None:
        for i in range(n):
            if i % points_per_line == 0:
                line_start[i] = True
        return line_start

    gaps = np.diff(flagged_idx)
    if len(gaps) == 0:
        return line_start
    median_gap = np.median(gaps)
    threshold = gap_factor * median_gap
    for i, g in enumerate(gaps, start=1):
        if g > threshold:
            line_start[i] = True
    return line_start


def generate_clocks(lines_csv, meas_csv, fs, points_per_line=None,
                     gap_factor=1.3, settle_delay_samples=0):
    """
    Build the full sample-indexed clock table.

    fs                    : DAC/sample rate in Hz (galvo waveform update rate)
    points_per_line       : force a fixed pixel count per line (e.g. 16 for
                             a 16x16 raster). Leave None to auto-detect from
                             sample-index gaps (works for nonuniform/
                             trapezoidal spacing too).
    gap_factor            : sensitivity of auto line-boundary detection.
    settle_delay_samples  : shift pixel_clk this many samples LATER than the
                             raw "arrived at target" flag, to allow galvo
                             mechanical settling before the pixel is
                             considered valid. 0 = use the flag as-is.

    Returns
    -------
    full_df   : one row per DAC sample (the master clock domain)
    edges_df  : one row per pixel (compact table, e.g. for TCSPC marker
                programming or quick lookup)
    """
    lines, meas = load_pair(lines_csv, meas_csv)
    n_total = len(lines)

    raw_valid = ((meas["a"] != 0) | (meas["b"] != 0)).to_numpy()
    flagged_idx = np.where(raw_valid)[0]

    if settle_delay_samples:
        shifted = np.clip(flagged_idx + settle_delay_samples, 0, n_total - 1)
        pixel_valid = np.zeros(n_total, dtype=bool)
        pixel_valid[shifted] = True
        flagged_idx = shifted
    else:
        pixel_valid = raw_valid

    n_pixels = len(flagged_idx)
    if n_pixels == 0:
        raise ValueError("No measurement points found (all rows are 0,0) in " + meas_csv)

    line_start_flags = detect_line_boundaries(flagged_idx, points_per_line, gap_factor)

    # Safety net: if auto-detection collapsed to a single line (e.g. gaps too
    # uniform for the chosen gap_factor) and the pixel count is a perfect
    # square, fall back to sqrt(n_pixels) pixels/line -- the natural
    # assumption for an N x N raster.
    if points_per_line is None and line_start_flags.sum() <= 1:
        root = int(round(n_pixels ** 0.5))
        if root * root == n_pixels:
            line_start_flags = detect_line_boundaries(flagged_idx, points_per_line=root)

    pixel_index_full = np.full(n_total, -1, dtype=int)
    line_index_full = np.full(n_total, -1, dtype=int)
    line_clk_full = np.zeros(n_total, dtype=bool)

    line_idx = -1
    for k, n in enumerate(flagged_idx):
        if line_start_flags[k]:
            line_idx += 1
            line_clk_full[n] = True
        pixel_index_full[n] = k
        line_index_full[n] = line_idx

    frame_clk_full = np.zeros(n_total, dtype=bool)
    frame_clk_full[0] = True
    frame_clk_full[-1] = True

    # pixel_valid_gate: high from this pixel's sample until just before the
    # next pixel's sample -- a "dwell window" rather than a single edge.
    # Useful for binning photon macro-times into the correct pixel.
    pixel_valid_gate = np.zeros(n_total, dtype=int)
    for k in range(n_pixels):
        start = flagged_idx[k]
        end = flagged_idx[k + 1] if k + 1 < n_pixels else n_total
        pixel_valid_gate[start:end] = k

    t = np.arange(n_total) / fs

    full_df = pd.DataFrame({
        "n": np.arange(n_total),
        "t_s": t,
        "x": lines["x"].to_numpy(),
        "y": lines["y"].to_numpy(),
        "pixel_clk": pixel_valid.astype(int),
        "pixel_index": pixel_index_full,
        "line_clk": line_clk_full.astype(int),
        "line_index": line_index_full,
        "frame_clk": frame_clk_full.astype(int),
        "pixel_valid_gate": pixel_valid_gate,
    })

    edges_df = pd.DataFrame({
        "pixel_index": np.arange(n_pixels),
        "sample_n": flagged_idx,
        "t_s": flagged_idx / fs,
        "x": lines["x"].to_numpy()[flagged_idx],
        "y": lines["y"].to_numpy()[flagged_idx],
        "line_index": line_index_full[flagged_idx],
        "is_line_start": line_start_flags.astype(int),
    })

    return full_df, edges_df


def write_summary(path, lines_csv, meas_csv, fs, full_df, edges_df):
    n_total = len(full_df)
    n_pixels = len(edges_df)
    n_lines = edges_df["line_index"].nunique()
    frame_duration = n_total / fs

    pixel_gaps = np.diff(edges_df["sample_n"].to_numpy())
    intra_line_gaps = pixel_gaps[edges_df["is_line_start"].to_numpy()[1:] == 0]
    line_break_gaps = pixel_gaps[edges_df["is_line_start"].to_numpy()[1:] == 1]

    lines_text = [
        f"Source lines file : {lines_csv}",
        f"Source meas file  : {meas_csv}",
        f"Sample rate (fs)  : {fs} Hz  (Ts = {1000/fs:.4f} ms)",
        f"Total samples     : {n_total}",
        f"Frame duration    : {frame_duration:.4f} s  ({1/frame_duration:.4f} Hz frame rate)",
        f"Total pixels      : {n_pixels}",
        f"Total lines       : {n_lines}  ({n_pixels/n_lines:.2f} pixels/line avg)",
        "",
        "Intra-line pixel spacing (samples):",
        f"  min={intra_line_gaps.min() if len(intra_line_gaps) else 'n/a'}  "
        f"max={intra_line_gaps.max() if len(intra_line_gaps) else 'n/a'}  "
        f"mean={intra_line_gaps.mean() if len(intra_line_gaps) else float('nan'):.2f}",
        f"  -> dwell time min={intra_line_gaps.min()/fs*1e3 if len(intra_line_gaps) else 0:.3f} ms, "
        f"max={intra_line_gaps.max()/fs*1e3 if len(intra_line_gaps) else 0:.3f} ms",
        "",
        "Line-turnaround (flyback) gap (samples):",
        f"  min={line_break_gaps.min() if len(line_break_gaps) else 'n/a'}  "
        f"max={line_break_gaps.max() if len(line_break_gaps) else 'n/a'}  "
        f"mean={line_break_gaps.mean() if len(line_break_gaps) else float('nan'):.2f}",
        f"  -> flyback time min={line_break_gaps.min()/fs*1e3 if len(line_break_gaps) else 0:.3f} ms, "
        f"max={line_break_gaps.max()/fs*1e3 if len(line_break_gaps) else 0:.3f} ms",
    ]
    text = "\n".join(lines_text)
    with open(path, "w") as f:
        f.write(text + "\n")
    print(text)


def plot_clocks(full_df, out_png, max_samples=1200):
    """Plot the first max_samples samples of pixel/line/frame clocks for a
    visual sanity check (full frame is usually too dense to read).
    Pass max_samples=None to plot the entire frame."""
    df = full_df if max_samples is None else full_df.iloc[:max_samples]
    fig, axes = plt.subplots(4, 1, figsize=(16, 9), sharex=True)

    axes[0].plot(df["t_s"], df["x"], label="X", linewidth=0.8)
    axes[0].plot(df["t_s"], df["y"], label="Y", linewidth=0.8)
    axes[0].set_ylabel("Galvo\nvoltage")
    axes[0].legend(loc="upper right", fontsize=8)
    axes[0].grid(True, alpha=0.3)

    axes[1].step(df["t_s"], df["pixel_clk"], where="post", color="tab:red")
    axes[1].set_ylabel("Pixel\nclock")
    axes[1].set_ylim(-0.1, 1.1)
    axes[1].grid(True, alpha=0.3)

    axes[2].step(df["t_s"], df["line_clk"], where="post", color="tab:green")
    axes[2].set_ylabel("Line\nclock")
    axes[2].set_ylim(-0.1, 1.1)
    axes[2].grid(True, alpha=0.3)

    axes[3].step(df["t_s"], df["frame_clk"], where="post", color="tab:purple")
    axes[3].set_ylabel("Frame\nclock")
    axes[3].set_ylim(-0.1, 1.1)
    axes[3].set_xlabel("Time (s)")
    axes[3].grid(True, alpha=0.3)

    label = "full frame" if max_samples is None else f"first {len(df)} samples"
    fig.suptitle(f"Sync clocks ({label}, {len(df)} of {len(full_df)} samples)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close(fig)


def run(lines_csv, meas_csv, fs, out_dir, points_per_line=None,
        gap_factor=1.5, settle_delay_samples=0):
    os.makedirs(out_dir, exist_ok=True)

    full_df, edges_df = generate_clocks(
        lines_csv, meas_csv, fs,
        points_per_line=points_per_line,
        gap_factor=gap_factor,
        settle_delay_samples=settle_delay_samples,
    )

    full_path = os.path.join(out_dir, "sync_clocks_full.csv")
    edges_path = os.path.join(out_dir, "sync_pixel_edges.csv")
    png_path = os.path.join(out_dir, "sync_clocks.png")
    png_full_path = os.path.join(out_dir, "sync_clocks_full.png")
    summary_path = os.path.join(out_dir, "summary.txt")

    full_df.to_csv(full_path, index=False)
    edges_df.to_csv(edges_path, index=False)
    plot_clocks(full_df, png_path)
    plot_clocks(full_df, png_full_path, max_samples=None)
    write_summary(summary_path, lines_csv, meas_csv, fs, full_df, edges_df)

    print(f"\nSaved:\n  {full_path}\n  {edges_path}\n  {png_path}\n  {png_full_path}\n  {summary_path}")
    return full_df, edges_df


def main():
    parser = argparse.ArgumentParser(description="Generate pixel/line/frame sync clocks from scan CSVs.")
    parser.add_argument("--lines", help="Path to *_lines.csv (x,y galvo voltage samples)")
    parser.add_argument("--meas", help="Path to *_meas_pts.csv (pixel flag per sample)")
    parser.add_argument("--fs", type=float, default=1000.0, help="Sample rate in Hz (default 1000)")
    parser.add_argument("--out", default="sync_out", help="Output directory")
    parser.add_argument("--points-per-line", type=int, default=None,
                         help="Force fixed pixels/line (e.g. 16). Default: auto-detect from gaps.")
    parser.add_argument("--gap-factor", type=float, default=1.3,
                         help="Auto line-break sensitivity (gap > factor * median intra-line gap)")
    parser.add_argument("--settle-delay-samples", type=int, default=0,
                         help="Delay pixel_clk this many samples for galvo settling")
    args = parser.parse_args()

    if args.lines and args.meas:
        run(args.lines, args.meas, args.fs, args.out,
            points_per_line=args.points_per_line,
            gap_factor=args.gap_factor,
            settle_delay_samples=args.settle_delay_samples)
        return

    # No args given: process both known 16x16 datasets as a convenience default.
    print("No --lines/--meas given, processing both default 16x16 datasets at fs=1000 Hz...\n")
    run(r"16x16_root_pair/root_lines.csv",
        r"16x16_root_pair/root_meas_pts.csv",
        fs=1000.0, out_dir="sync_out_root")
    print()
    run(r"16x16_square_Pattern/16x16_lines.csv",
        r"16x16_square_Pattern/meas_pts_16x.csv",
        fs=1000.0, out_dir="sync_out_square")


if __name__ == "__main__":
    main()
