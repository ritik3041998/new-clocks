"""
clocks_from_equal_spacing.py

Builds pixel/line/frame sync clocks from the output of
equal_spacing_from_trajectory.py (equal_spaced_pixels.csv), and renders
both the full clock-timeline plot and a trajectory-with-clock-events plot.

Usage
-----
    python clocks_from_equal_spacing.py \
        --lines <lines.csv> --pixels <equal_spaced_pixels.csv> \
        --fs 1000 --out <out_dir> --title "root_pair"
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def build_clock_table(lines_csv, pixels_csv, fs):
    lines = pd.read_csv(lines_csv, header=None, names=["x", "y"])
    pixels = pd.read_csv(pixels_csv)
    n_total = len(lines)

    sample_n = pixels["sample_n"].to_numpy()
    line_index = pixels["line_index"].to_numpy()
    n_lines = line_index.max() + 1

    pixel_clk = np.zeros(n_total, dtype=int)
    pixel_index_full = np.full(n_total, -1, dtype=int)
    line_index_full = np.full(n_total, -1, dtype=int)
    line_clk = np.zeros(n_total, dtype=int)

    pixel_clk[sample_n] = 1
    pixel_index_full[sample_n] = np.arange(len(sample_n))
    line_index_full[sample_n] = line_index

    # line_clk fires on the first pixel of each line
    first_of_line = pixels.groupby("line_index")["sample_n"].min().to_numpy()
    line_clk[first_of_line] = 1

    # frame_clk brackets the actual pixel data -- fires at the first and
    # last detected pixel_clk sample, not at the raw buffer edges (which
    # may include lead-in/settle samples before the first real pixel and
    # trailing samples after the last one).
    frame_clk = np.zeros(n_total, dtype=int)
    frame_clk[sample_n.min()] = 1
    frame_clk[sample_n.max()] = 1

    t = np.arange(n_total) / fs

    full_df = pd.DataFrame({
        "n": np.arange(n_total),
        "t_s": t,
        "x": lines["x"].to_numpy(),
        "y": lines["y"].to_numpy(),
        "pixel_clk": pixel_clk,
        "pixel_index": pixel_index_full,
        "line_clk": line_clk,
        "line_index": line_index_full,
        "frame_clk": frame_clk,
    })
    return full_df, pixels, n_lines


def write_summary(path, full_df, pixels, fs, n_lines):
    n_total = len(full_df)
    n_pixels = len(pixels)
    frame_duration = n_total / fs

    gaps = np.diff(pixels["sample_n"].to_numpy())
    is_first = (pixels["line_index"].diff().fillna(1) != 0).to_numpy()[1:]
    intra = gaps[~is_first]
    inter = gaps[is_first]

    lines_out = [
        f"Sample rate (fs) : {fs} Hz  (Ts = {1000/fs:.4f} ms)",
        f"Total samples    : {n_total}",
        f"Frame duration   : {frame_duration:.4f} s  ({1/frame_duration:.4f} Hz frame rate)",
        f"Total pixels     : {n_pixels}",
        f"Total lines      : {n_lines}  ({n_pixels/n_lines:.2f} pixels/line)",
        "",
        "Intra-line pixel spacing (samples):",
        f"  min={intra.min() if len(intra) else 'n/a'}  max={intra.max() if len(intra) else 'n/a'}  "
        f"mean={intra.mean() if len(intra) else float('nan'):.2f}",
        "",
        "Line-turnaround gap (samples):",
        f"  min={inter.min() if len(inter) else 'n/a'}  max={inter.max() if len(inter) else 'n/a'}  "
        f"mean={inter.mean() if len(inter) else float('nan'):.2f}",
    ]
    text = "\n".join(lines_out)
    with open(path, "w") as f:
        f.write(text + "\n")
    print(text)


def build_synthetic_clock_table(lines_csv, pixels_csv, fs):
    """
    Build the pixel/line/frame clock timeline WITHOUT moving any of the
    pixel dot positions in pixels_csv. The dot positions (x,y) are already
    correctly equally-spaced per line and are left untouched (the
    trajectory plot built from pixels_csv is unaffected by this function).

    What this DOES fix: the *timing* of those same pixels. Each line's own
    arc length (already correct, from the unchanged x,y) is used to derive
    a per-line clock duration at a constant scan velocity -- so a shorter
    line's 16 (already-correct) points get a proportionally shorter pixel-
    clock burst, exactly mirroring how close together they already are in
    space. Line order and inter-line gaps come from the original data;
    only the *spacing in time* of the clock pulses within each burst
    changes to track distance instead of being a fixed constant.
    """
    lines = pd.read_csv(lines_csv, header=None, names=["x", "y"])
    pixels = pd.read_csv(pixels_csv)
    n_lines = pixels["line_index"].max() + 1
    points_per_line = len(pixels) // n_lines

    # Reference: line 0's own already-detected duration and arc length
    # (from the unchanged data) set the scan velocity (voltage units/sample).
    def line_arc_length(seg):
        return np.sqrt(np.diff(seg["x"]) ** 2 + np.diff(seg["y"]) ** 2).sum()

    seg0 = pixels[pixels.line_index == 0]
    dur0 = seg0["sample_n"].max() - seg0["sample_n"].min()
    len0 = line_arc_length(seg0)
    velocity = len0 / dur0  # voltage units per sample

    # Typical inter-line gap in the original (unchanged) data -- reused
    # as-is, since gap/flyback timing isn't part of this fix.
    starts_orig = pixels.groupby("line_index")["sample_n"].min().to_numpy()
    ends_orig = pixels.groupby("line_index")["sample_n"].max().to_numpy()
    gaps = starts_orig[1:] - ends_orig[:-1]
    base_gap = int(round(np.median(gaps))) if len(gaps) else 44

    lead_in = int(starts_orig[0])  # keep the same start offset as the original

    pixel_events = []  # (sample_n, line_index)
    line_starts = []
    cursor = lead_in
    for L in range(n_lines):
        seg = pixels[pixels.line_index == L]
        arc = line_arc_length(seg)
        duration = max(points_per_line - 1, int(round(arc / velocity)))
        line_starts.append(cursor)
        times = np.linspace(cursor, cursor + duration, points_per_line)
        times = np.round(times).astype(int)
        for t in times:
            pixel_events.append((t, L))
        cursor = cursor + duration + base_gap

    trail = lead_in  # mirror the same amount of trailing buffer as lead-in
    n_total_synthetic = cursor + trail
    n_total = max(n_total_synthetic, len(lines))

    pixel_clk = np.zeros(n_total, dtype=int)
    pixel_index_full = np.full(n_total, -1, dtype=int)
    line_index_full = np.full(n_total, -1, dtype=int)
    line_clk = np.zeros(n_total, dtype=int)

    for i, (t, L) in enumerate(pixel_events):
        pixel_clk[t] = 1
        pixel_index_full[t] = i
        line_index_full[t] = L
    for L, s in enumerate(line_starts):
        line_clk[s] = 1

    # frame_clk brackets the actual pixel events -- fires at the first and
    # last detected pixel_clk sample, not at the raw buffer edges.
    pixel_sample_ns = [t for t, _ in pixel_events]
    frame_clk = np.zeros(n_total, dtype=int)
    frame_clk[min(pixel_sample_ns)] = 1
    frame_clk[max(pixel_sample_ns)] = 1

    # X,Y kept from the original trajectory for visual context (padded/
    # truncated to the synthetic length -- purely a background reference,
    # the clock signals above are the point of this table).
    x_full = lines["x"].to_numpy()
    y_full = lines["y"].to_numpy()
    if len(x_full) < n_total:
        x_full = np.pad(x_full, (0, n_total - len(x_full)), mode="edge")
        y_full = np.pad(y_full, (0, n_total - len(y_full)), mode="edge")

    t = np.arange(n_total) / fs
    full_df = pd.DataFrame({
        "n": np.arange(n_total),
        "t_s": t,
        "x": x_full[:n_total],
        "y": y_full[:n_total],
        "pixel_clk": pixel_clk,
        "pixel_index": pixel_index_full,
        "line_clk": line_clk,
        "line_index": line_index_full,
        "frame_clk": frame_clk,
    })
    return full_df


def plot_clock_timeline(full_df, out_png, title):
    fig, axes = plt.subplots(3, 1, figsize=(16, 7), sharex=True)

    axes[0].step(full_df["t_s"], full_df["pixel_clk"], where="post", color="tab:red")
    axes[0].set_ylabel("Pixel\nclock")
    axes[0].set_ylim(-0.1, 1.1)
    axes[0].grid(True, alpha=0.3)

    axes[1].step(full_df["t_s"], full_df["line_clk"], where="post", color="tab:green")
    axes[1].set_ylabel("Line\nclock")
    axes[1].set_ylim(-0.1, 1.1)
    axes[1].grid(True, alpha=0.3)

    axes[2].step(full_df["t_s"], full_df["frame_clk"], where="post", color="tab:purple")
    axes[2].set_ylabel("Frame\nclock")
    axes[2].set_ylim(-0.1, 1.1)
    axes[2].set_xlabel("Time (s)")
    axes[2].grid(True, alpha=0.3)

    fig.suptitle(f"{title}: sync clocks from equal-spaced pixels ({len(full_df)} samples)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_trajectory_with_events(full_df, out_png, title, dark=True):
    fig, ax = plt.subplots(figsize=(7, 7))
    if dark:
        fig.patch.set_facecolor("black")
        ax.set_facecolor("black")
        traj_color, text_color, pixel_color, line_color = "white", "white", "red", "cyan"
    else:
        traj_color, text_color, pixel_color, line_color = "tab:blue", "black", "red", "tab:green"

    ax.plot(full_df["x"], full_df["y"], color=traj_color, linewidth=0.8, zorder=1)
    pix = full_df[full_df["pixel_clk"] == 1]
    ax.scatter(pix["x"], pix["y"], s=16, color=pixel_color, zorder=2, label="Pixel clock")
    linept = full_df[full_df["line_clk"] == 1]
    ax.scatter(linept["x"], linept["y"], s=40, facecolors="none", edgecolors=line_color,
               linewidths=1.4, zorder=3, label="Line clock (first pixel)")

    ax.set_title(f"{title}: Scan Trajectory with Pixel/Line-Clock Events", color=text_color)
    ax.set_xlabel("X Voltage", color=text_color)
    ax.set_ylabel("Y Voltage", color=text_color)
    ax.tick_params(colors=text_color)
    for spine in ax.spines.values():
        spine.set_color(text_color)
    legend = ax.legend(fontsize=8, loc="upper right")
    if dark:
        legend.get_frame().set_facecolor("black")
        for t in legend.get_texts():
            t.set_color("white")
    ax.axis("equal")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lines", required=True)
    parser.add_argument("--pixels", required=True, help="equal_spaced_pixels.csv from equal_spacing_from_trajectory.py")
    parser.add_argument("--fs", type=float, default=1000.0)
    parser.add_argument("--out", default="clocks_out")
    parser.add_argument("--title", default="pattern")
    parser.add_argument("--light", action="store_true", help="Use light theme instead of dark for the trajectory plot")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    full_df, pixels, n_lines = build_clock_table(args.lines, args.pixels, args.fs)

    full_csv = os.path.join(args.out, "sync_clocks_full.csv")
    timeline_png = os.path.join(args.out, "sync_clocks_timeline.png")
    traj_png = os.path.join(args.out, "sync_clocks_trajectory.png")
    summary_txt = os.path.join(args.out, "summary.txt")

    full_df.to_csv(full_csv, index=False)
    plot_clock_timeline(full_df, timeline_png, args.title)
    plot_trajectory_with_events(full_df, traj_png, args.title, dark=not args.light)
    write_summary(summary_txt, full_df, pixels, args.fs, n_lines)

    print(f"\nSaved:\n  {full_csv}\n  {timeline_png}\n  {traj_png}\n  {summary_txt}")


if __name__ == "__main__":
    main()
