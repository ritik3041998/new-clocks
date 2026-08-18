"""
equal_spacing_from_trajectory.py

Places N equally arc-length-spaced pixel markers on each line of a raster
scan trajectory, using ONLY the x,y voltage trace (no separate meas_pts.csv
required).

Method (v2 -- curvature-based straight-run detection)
-------------------------------------------------------
The first version of this script guessed each line's active window by
trimming a fixed fraction off both ends of a turn-to-turn segment. That
fraction doesn't scale correctly across lines of very different physical
length (a trapezoidal/corrected pattern has a different line length on
every line), so short lines ended up with points crammed into a fraction
of the visible line, leaving the rest empty.

This version instead finds the actual STRAIGHT portion of each line
directly, using path curvature:
  1. Compute the local turning angle at every sample (how much the
     direction vector rotates from one sample to the next).
  2. Smooth it, then classify samples as "straight" (low curvature) or
     "turning" (high curvature, i.e. inside a U-turn/flyback).
  3. Take the contiguous straight runs, keep the n_lines longest ones
     (sorted back into chronological/scan order) -- these ARE the lines,
     found directly from geometry, no fixed trim fraction needed.
  4. Within each straight run's full extent (start to end, no guessed
     margin), place points_per_line points at equal arc-length spacing --
     the same guarantee as before, but now covering the FULL line just
     like the direct/uniform square pattern does.

Usage
-----
    python equal_spacing_from_trajectory.py --lines <lines.csv> \
        --n-lines 16 --points-per-line 16 --out out_dir \
        [--compare-meas <meas_pts.csv>]   # optional ground-truth check
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def compute_curvature(x, y, smooth_window=5):
    dx = np.diff(x)
    dy = np.diff(y)
    ang = np.arctan2(dy, dx)
    dang = np.diff(ang)
    # wrap to [-pi, pi]
    dang = (dang + np.pi) % (2 * np.pi) - np.pi
    kappa = np.abs(dang)
    kappa = np.concatenate([[kappa[0]], kappa, [kappa[-1]]])  # pad back to len(x)
    if smooth_window > 1:
        kernel = np.ones(smooth_window) / smooth_window
        kappa = np.convolve(kappa, kernel, mode="same")
    return kappa


def find_straight_runs(kappa, n_lines, min_run_frac=0.3, straight_percentile=92):
    """
    Classify samples as straight/turning by thresholding curvature, find
    contiguous straight runs, and keep the n_lines longest ones (returned
    in chronological order = scan order).
    """
    n = len(kappa)
    # Threshold: turns have much higher curvature than straight runs, but
    # they are SHORT arcs -- most samples in the whole trajectory are on a
    # straight run. A high percentile (~90-95) correctly keeps almost all
    # straight samples while excluding the turn arcs. Too low a percentile
    # (e.g. 70) starts cutting into the straight runs themselves near their
    # ends, fragmenting them and shrinking the usable line length.
    threshold = np.percentile(kappa, straight_percentile)
    # Use <= so a perfectly straight trajectory (curvature exactly 0 on
    # straight runs, e.g. the direct square scan) still classifies as
    # straight even when the percentile threshold itself lands on 0.
    is_straight = kappa <= threshold

    # find contiguous True runs
    runs = []
    start = None
    for i, v in enumerate(is_straight):
        if v and start is None:
            start = i
        elif not v and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, n - 1))

    min_len = n / n_lines * min_run_frac
    runs = [r for r in runs if (r[1] - r[0]) >= min_len]

    if len(runs) < n_lines:
        raise ValueError(
            f"Found only {len(runs)} straight runs >= {min_len:.0f} samples, "
            f"need {n_lines}. Try lowering --min-run-frac."
        )

    # keep the n_lines longest, then re-sort chronologically
    runs_by_len = sorted(runs, key=lambda r: r[1] - r[0], reverse=True)[:n_lines]
    runs_sorted = sorted(runs_by_len, key=lambda r: r[0])
    return runs_sorted


def place_points_on_segment(x, y, start, end, points_per_line, edge_trim_samples=1):
    """
    Within sample range [start, end] (the full straight run), trim a small
    fixed number of samples off each end (numerical safety only -- not a
    fraction of line length), then place points_per_line points at equal
    ARC-LENGTH spacing across the (nearly) full remaining span.
    """
    lo = start + edge_trim_samples
    hi = end - edge_trim_samples
    if hi <= lo:
        lo, hi = start, end

    xs = x[lo:hi + 1]
    ys = y[lo:hi + 1]
    d = np.sqrt(np.diff(xs) ** 2 + np.diff(ys) ** 2)
    cum = np.concatenate([[0], np.cumsum(d)])
    total = cum[-1]

    target_s = np.linspace(0, total, points_per_line)
    exact_x = np.interp(target_s, cum, xs)
    exact_y = np.interp(target_s, cum, ys)
    sample_indices = np.interp(target_s, cum, np.arange(lo, hi + 1))
    sample_indices = np.round(sample_indices).astype(int)
    sample_indices = np.clip(sample_indices, 0, len(x) - 1)
    return sample_indices, exact_x, exact_y


def arc_length(x, y, start, end):
    xs = x[start:end + 1]
    ys = y[start:end + 1]
    return np.sqrt(np.diff(xs) ** 2 + np.diff(ys) ** 2).sum()


def rescale_segments_by_velocity(x, y, segments):
    """
    The curvature-threshold detector finds each line's START reliably, but
    can over-extend the END into a low-curvature tail of the turn, which
    hides the fact that shorter (narrower) lines should take FEWER samples
    -- distance and duration should shrink together, at a roughly constant
    scan velocity (voltage units per sample), the same way the ground-truth
    ...meas_pts.csv data behaves.

    Fix: compute each segment's arc length, take the LONGEST segment as the
    reference (least ambiguous -- a wide line's true end is easiest to find),
    derive velocity = reference arc length / reference duration, then trim
    every other segment's END so its sample count matches its own arc
    length at that same velocity. Shorter lines end up with visibly shorter,
    denser sample windows, exactly mirroring how their points are already
    packed closer together in the trajectory plot.
    """
    lengths = [arc_length(x, y, s, e) for s, e in segments]
    durations = [e - s for s, e in segments]
    ref = int(np.argmax(lengths))
    velocity = lengths[ref] / durations[ref]  # voltage units per sample

    rescaled = []
    for (s, e), L in zip(segments, lengths):
        target_dur = max(1, int(round(L / velocity)))
        new_e = min(e, s + target_dur)
        rescaled.append((s, new_e))
    return rescaled


def generate(lines_csv, n_lines, points_per_line, min_run_frac=0.3, edge_trim_samples=1,
             smooth_window=5, straight_percentile=92, scale_by_velocity=False):
    lines = pd.read_csv(lines_csv, header=None, names=["x", "y"])
    x = lines["x"].to_numpy()
    y = lines["y"].to_numpy()

    kappa = compute_curvature(x, y, smooth_window)
    segments = find_straight_runs(kappa, n_lines, min_run_frac, straight_percentile)
    if scale_by_velocity:
        segments = rescale_segments_by_velocity(x, y, segments)

    all_pixel_idx, all_exact_x, all_exact_y, line_of_pixel = [], [], [], []
    for L, (start, end) in enumerate(segments):
        idxs, exact_x, exact_y = place_points_on_segment(x, y, start, end, points_per_line, edge_trim_samples)
        all_pixel_idx.extend(idxs.tolist())
        all_exact_x.extend(exact_x.tolist())
        all_exact_y.extend(exact_y.tolist())
        line_of_pixel.extend([L] * len(idxs))

    edges_df = pd.DataFrame({
        "pixel_index": np.arange(len(all_pixel_idx)),
        "sample_n": all_pixel_idx,
        "x": all_exact_x,
        "y": all_exact_y,
        "x_nearest_sample": x[all_pixel_idx],
        "y_nearest_sample": y[all_pixel_idx],
        "line_index": line_of_pixel,
    })
    return lines, edges_df, segments


def compare_to_ground_truth(edges_df, meas_csv, lines_csv):
    lines = pd.read_csv(lines_csv, header=None, names=["x", "y"])
    meas = pd.read_csv(meas_csv, header=None, names=["a", "b"])
    gt_idx = meas.index[(meas["a"] != 0) | (meas["b"] != 0)].to_numpy()
    gt_pts = lines.loc[gt_idx].reset_index(drop=True)

    n = min(len(gt_pts), len(edges_df))
    dx = edges_df["x"].to_numpy()[:n] - gt_pts["x"].to_numpy()[:n]
    dy = edges_df["y"].to_numpy()[:n] - gt_pts["y"].to_numpy()[:n]
    err = np.sqrt(dx**2 + dy**2)
    print(f"\nValidation against ground-truth meas file ({meas_csv}):")
    print(f"  matched points : {n}")
    print(f"  position error : mean={err.mean():.6f}  max={err.max():.6f}  "
          f"(voltage units -- compare to typical line pitch ~0.1-0.26)")


def spacing_report(edges_df, n_lines):
    print(f"\nPer-line spacing uniformity ({n_lines} lines):")
    print(f"{'line':>4} {'n_pts':>6} {'mean spacing':>13} {'std dev':>10} {'%variation':>11}")
    for L in range(n_lines):
        seg = edges_df[edges_df.line_index == L]
        d = np.sqrt(np.diff(seg["x"]) ** 2 + np.diff(seg["y"]) ** 2)
        pct = 100 * d.std() / d.mean() if d.mean() > 0 else 0
        print(f"{L:>4} {len(seg):>6} {d.mean():>13.6f} {d.std():>10.8f} {pct:>10.4f}%")


def plot_result(lines, edges_df, out_png, segments=None):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(lines["x"], lines["y"], linewidth=0.7, color="tab:blue", zorder=1, label="Trajectory")
    ax.scatter(edges_df["x"], edges_df["y"], s=14, color="red", zorder=2, label="Equally-spaced pixels")
    ax.set_title(f"Equal arc-length spacing per line ({len(edges_df)} points)")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.axis("equal")
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Place equally arc-length-spaced pixels per line, from trajectory only.")
    parser.add_argument("--lines", required=True, help="Path to lines.csv (x,y galvo voltage samples)")
    parser.add_argument("--n-lines", type=int, default=16)
    parser.add_argument("--points-per-line", type=int, default=16)
    parser.add_argument("--min-run-frac", type=float, default=0.3,
                         help="Minimum straight-run length as a fraction of the expected per-line sample count")
    parser.add_argument("--edge-trim-samples", type=int, default=1,
                         help="Fixed number of samples trimmed off each straight-run end (numerical safety only)")
    parser.add_argument("--smooth-window", type=int, default=5,
                         help="Curvature smoothing window (samples)")
    parser.add_argument("--straight-percentile", type=float, default=92,
                         help="Percentile of curvature below which a sample counts as 'straight' (higher = more permissive)")
    parser.add_argument("--out", default="equal_spacing_out")
    parser.add_argument("--compare-meas", default=None,
                         help="Optional ground-truth meas_pts.csv to validate against")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    lines, edges_df, segments = generate(
        args.lines, args.n_lines, args.points_per_line,
        args.min_run_frac, args.edge_trim_samples, args.smooth_window,
        args.straight_percentile,
    )

    edges_path = os.path.join(args.out, "equal_spaced_pixels.csv")
    png_path = os.path.join(args.out, "equal_spacing_check.png")
    edges_df.to_csv(edges_path, index=False)
    plot_result(lines, edges_df, png_path, segments)

    print(f"Detected {len(segments)} straight-run line segments, placed "
          f"{args.points_per_line} points each ({len(edges_df)} total pixels)")
    print(f"Segment sample lengths: {[e - s for s, e in segments]}")
    print(f"Saved:\n  {edges_path}\n  {png_path}")

    spacing_report(edges_df, args.n_lines)

    if args.compare_meas:
        compare_to_ground_truth(edges_df, args.compare_meas, args.lines)


if __name__ == "__main__":
    main()
