# Galvo Scan Sync-Clock Generation — Workflow & Reference

This documents the complete pipeline built in this project: turning a galvo
scan trajectory (X/Y voltage samples) into **pixel / line / frame sync
clocks** for a TCSPC system, for both the direct **square** scan and the
geometrically-**corrected trapezoidal** scan (the "root_pair" pattern used
to counteract keystone/pincushion distortion from a tilted NLOS wall).

---

## 1. The data model

Every pattern is described by two parallel, sample-aligned CSVs (no header,
row *n* in one corresponds to row *n* in the other):

| File | Columns | Meaning |
|---|---|---|
| `*_lines.csv` | `x, y` | The galvo voltage commanded at DAC sample *n*. This is the full trajectory — thousands of samples tracing the serpentine/raster path. |
| `*_meas_pts.csv` / `meas_pts_16x.csv` | `a, b` | A **flag**, not a coordinate: `(0,0)` = not a pixel, `(5,5)` = "a measurement happens at this sample." Exactly 256 rows are flagged for a 16×16 pattern. |

Both patterns in this project have `8349` total samples and `256` flagged
pixel samples (16 lines × 16 pixels).

**Why two files instead of one:** the trajectory (`lines.csv`) is a
continuous curve — it also includes turnarounds, settle time, and flyback.
The flag file marks only the handful of samples where the galvo is actually
sitting *on* a grid point, ready for a measurement. Clock generation is
built entirely on that flag file when it's available — never re-derived
from scratch — because it's the authoritative source.

---

## 2. Two different situations, two different pipelines

| Situation | Script | Trust level |
|---|---|---|
| You **have** the real `*_meas_pts.csv` (this project's given data) | [`generate_sync_clocks.py`](generate_sync_clocks.py) | **Exact.** Reads the flags directly — no reconstruction, no approximation. |
| You only have the raw trajectory (`*_lines.csv`), e.g. for a **brand-new pattern** with no flag file yet | [`equal_spacing_from_trajectory.py`](equal_spacing_from_trajectory.py) → [`clocks_from_equal_spacing.py`](clocks_from_equal_spacing.py) | **Approximate but principled.** Reconstructs pixel positions and timing directly from trajectory geometry. |

Both pipelines produce the same three outputs: a `pixel_clk`, a `line_clk`,
and a `frame_clk`, all indexed by the same master sample counter `n` that
also indexes the DAC voltage samples — this is what keeps everything
deterministic and jitter-free when ported to real hardware (FPGA/DAQ).

---

## 3. Pipeline A — ground truth (`generate_sync_clocks.py`)

This is the simple, authoritative case, used whenever `*_meas_pts.csv`
exists.

```
python generate_sync_clocks.py --lines <lines.csv> --meas <meas_pts.csv> --fs 1000 --out <out_dir>
```

**Steps:**
1. Load `lines.csv` (voltage samples) and `meas_pts.csv` (flags) — same
   row count, sample-aligned.
2. `pixel_clk[n] = 1` wherever the flag is nonzero. That's it — no
   detection, no guessing.
3. **Line boundaries** are found from the *gaps* between consecutive
   flagged samples: within a line, gaps are small and consistent (e.g. 30
   samples); at a flyback/turnaround, the gap jumps sharply (e.g. 44+
   samples). A gap `> 1.3 × median gap` marks the first pixel of a new
   line → `line_clk[n] = 1` there.
4. `frame_clk[n] = 1` at `n = 0` and `n = N_total − 1`.
5. Every sample also gets a `pixel_valid_gate` — the *dwell window* from
   one pixel clock to the next, useful for binning photon macro-times into
   the right pixel rather than relying on instantaneous edges alone.
6. Outputs: `sync_clocks_full.csv` (per-sample table), `sync_pixel_edges.csv`
   (one row per pixel), a full-frame timeline PNG, and `summary.txt` with
   dwell/flyback statistics.

**Why this is exact for the trapezoidal case:** the flag file was generated
by the same process that computed the corrected voltages, so it already
encodes the true, physically-shrinking per-line pitch. Verified against the
original `scan_trajectory_dark.png` reference — pixel-to-voltage ratio
matched to within 0.3% on every one of the 16 lines.

---

## 4. Pipeline B — reconstruction from trajectory only

Used when there is **no** `meas_pts.csv` — e.g. designing sync clocks for a
new pattern before the measurement-point flags exist yet. Two scripts, two
separate jobs:

### 4a. Where are the pixels? — `equal_spacing_from_trajectory.py`

**Goal:** place `points_per_line` pixels on each of `n_lines` scan lines,
evenly spaced by *arc length* (physical distance along the path), reading
only `x, y` from `lines.csv`.

**Step 1 — find where each line actually is, using curvature, not a guess.**
A first version of this script tried to find line boundaries via Y-axis
extrema (peaks/valleys) and trim a *fixed percentage* off both ends as a
settle margin. That failed for the trapezoidal case: a fixed fraction
doesn't scale to lines of very different physical length, so short lines
got their 16 points crammed into a fraction of the visible line, leaving
the rest empty.

The fix computes the **local turning angle** at every sample:

```
kappa[n] = |angle between direction(n-1→n) and direction(n→n+1)|
```

Straight runs have `kappa ≈ 0`; the rounded U-turns/flybacks have high
`kappa`. Samples are classified "straight" if `kappa` is below the 92nd
percentile of the whole trajectory's curvature — high because turns are
short arcs, so the overwhelming majority of samples on any reasonable scan
*are* straight, and the threshold should barely touch them. The 16 longest
contiguous straight runs, re-sorted back into chronological order, **are**
the lines — found from geometry, not from an assumed shape.

**Step 2 — place points at equal arc length within each run.**
For each detected line, cumulative Euclidean distance is computed sample by
sample; `points_per_line` targets are placed at equal intervals of that
cumulative distance, then mapped back to `(x, y)` by linear interpolation
(exact, continuous — not snapped to the nearest existing sample). Result:
**0.0000% spacing variation** within every line, verified on both the
square and trapezoidal datasets, with full edge-to-edge coverage (no more
empty gaps).

This step is intentionally **only about where the dots are** — it does not
decide timing. `equal_spaced_pixels.csv` stores each pixel's `(x, y)`,
`line_index`, and the *nearest existing sample* (`sample_n`) it lands
closest to.

### 4b. When does each pixel happen? — `clocks_from_equal_spacing.py`

Turns the pixel table from 4a into a per-sample clock table, timeline PNG,
and trajectory-with-events PNG. Two different timing modes exist here,
because they answer different questions:

**Direct mode (`build_clock_table`)** — uses each pixel's `sample_n` as-is.
This is correct for the **square pattern**, where the trajectory really
does move at constant sample-rate-per-voltage everywhere, so `sample_n`
already reflects true timing (verified: matches ground truth to ~0.0045
average position error).

**Synthetic mode (`build_synthetic_clock_table`)** — used for the
**trapezoidal pattern**, and this is the part worth explaining carefully,
because it's the crux of "how do you get a trapezoidal *clock*, not just a
trapezoidal *shape*":

> The curvature detector in step 4a is very good at finding where each line
> **starts**, but can over-extend where it thinks a line **ends** (it
> tolerates a little residual curvature at the tail of a turn). The result:
> every detected line came out ≈452 samples long regardless of its true
> physical length — hiding the fact that a physically shorter (narrower)
> line should also take *less time* to scan, not just less space.

The fix does **not** touch the pixel positions from 4a at all (they're
already correct and are left completely alone). It only recomputes *when*
each of those same pixels fires, using this reasoning:

1. Take line 0 (the widest, least ambiguous line) as a **velocity
   reference**: `velocity = arc_length(line 0) / duration(line 0)` in
   voltage-units per sample — i.e. how fast the galvo moves along a
   straight run.
2. For every other line, its arc length is already known exactly (from the
   unchanged, already-correct point positions). Its **duration** is
   recomputed as `arc_length(line L) / velocity` — the same constant
   velocity, applied to that line's own (shorter) physical length. A
   narrower line finishes proportionally faster.
3. The 16 (still-unmoved) pixel positions of that line are laid out at
   **equal time intervals** across this newly-computed, shorter duration —
   `line_clk` fires at the start, `pixel_clk` fires at each of the 16
   evenly-spaced instants, and the next line begins after the same
   flyback gap used elsewhere in the data.
4. `frame_clk` brackets the whole synthetic sequence.

Net effect: the **timeline** now visibly narrows toward the later
(narrower) lines — burst duration shrinks from ~462 ms down to ~361 ms
across the 16 lines in this dataset — exactly mirroring how those same
pixels are already packed closer together in space, without moving a
single dot in the trajectory plot.

```
python equal_spacing_from_trajectory.py --lines <lines.csv> --n-lines 16 --points-per-line 16 --out <eq_dir>
python clocks_from_equal_spacing.py --lines <lines.csv> --pixels <eq_dir>/equal_spaced_pixels.csv --fs 1000 --out <clk_dir> --title "<name>"
```

---

## 5. How the trapezoidal *shape* becomes a trapezoidal *clock*, end to end

Putting §4a and §4b together, here is the full causal chain for the
corrected/trapezoidal pattern specifically:

```
Geometric correction (external, already given)
   → produces X,Y voltages that trace a trapezoid in galvo-voltage space
     so the projection lands as a uniform grid on the tilted wall
        │
        ▼
lines.csv (8349 samples, continuous trajectory)
        │
        ▼  curvature-based straight-run detection (equal_spacing_from_trajectory.py)
16 line segments, each with its OWN arc length
(shorter for narrower lines near the trapezoid's compressed edge)
        │
        ▼  equal arc-length placement within each segment
256 pixel positions, 0% spacing variation per line, full coverage
(equal_spaced_pixels.csv -- SPATIAL correctness, done)
        │
        ▼  velocity-based duration rescaling (clocks_from_equal_spacing.py)
16 line DURATIONS, each proportional to that line's own arc length
at one constant reference velocity
        │
        ▼  lay out pixel/line/frame events on a synthetic per-sample timeline
pixel_clk / line_clk / frame_clk
(TEMPORAL correctness: narrower lines -> shorter bursts -> higher pixel-clock frequency)
```

The **spatial** question ("where are the 256 points") and the **temporal**
question ("when does each one fire") are solved by two separate, composable
steps — this is deliberate, matching the general design principle from the
timing-design analysis earlier in this project: geometric/spatial
correction and photon-timing synchronization are different domains and
should never be conflated into one computation.

---

## 6. Verifying against ground truth

Whenever a real `meas_pts.csv` is available, `equal_spacing_from_trajectory.py --compare-meas <file>`
reports the average/max position error between the reconstruction and the
real flagged points — this project's trapezoidal reconstruction lands at
≈0.23–0.42 voltage-units average error (a real approximation, useful when
no ground truth exists), while the square-pattern reconstruction lands at
≈0.0045 (essentially exact, because the square scan's assumptions —
constant velocity, constant line length — are actually true for it).

**Rule of thumb used throughout this project:** if a real `meas_pts.csv`
exists for the pattern you're generating clocks for, always prefer Pipeline
A. Pipeline B exists specifically for the case where it doesn't.

---

## 7. Every parameter produced by clock generation

Both pipelines write a `sync_clocks_full.csv` — one row per DAC sample,
the master timeline everything else is derived from. This is the complete
column reference:

| Column | Type | Meaning |
|---|---|---|
| `n` | int | The master sample counter, `0 … N_total−1`. This is the single index that ties the DAC voltage, the pixel clock, the line clock, and the frame clock together — the reason the whole scheme is jitter-free. |
| `t_s` | float | Wall-clock time of sample `n`, i.e. `n / fs`. Purely a display/derived convenience — all real logic runs on `n`, not `t_s`, so changing `fs` never changes any other column. |
| `x`, `y` | float | The commanded galvo voltage at sample `n` (from `lines.csv`). Kept in the table for plotting/debugging context; not needed by the TCSPC. |
| `pixel_clk` | 0/1 | **1** on samples where a pixel measurement should be taken, else 0. This is the signal that goes to the TCSPC's pixel/marker input. |
| `pixel_index` | int (−1 if none) | Which of the 256 pixels this `pixel_clk` edge is, in scan order `0 … 255`. Lets you map a marker event straight back to `(line, column)` without recomputing anything. |
| `line_clk` | 0/1 | **1** on the sample that is the *first pixel of a new line* — one edge per line (16 per frame here). Goes to the TCSPC's line/marker input, or increments a line counter downstream. |
| `line_index` | int (−1 if none) | Which of the 16 lines this sample's pixel belongs to, `0 … 15`. |
| `frame_clk` | 0/1 | **1** exactly twice per frame: at the first pixel sample and the last pixel sample (see §8 below — deliberately *not* the raw buffer edges, which can include lead-in/settle dead time). |
| `pixel_valid_gate` *(ground-truth pipeline only)* | int | The **dwell window** index — constant at pixel index `k` from that pixel's `pixel_clk` sample up to (but not including) the next one. This is what you actually want for *binning* a photon: "which window does this photon's arrival time fall into," rather than only having instantaneous edges. |

`sync_pixel_edges.csv` / `equal_spaced_pixels.csv` are the compact,
one-row-per-pixel companions — `pixel_index, sample_n, t_s, x, y,
line_index, is_line_start` — useful for pre-loading a TCSPC's marker/ROI
table or for a host-side lookup without scanning the full per-sample array.

`summary.txt` reports, per pattern: total samples, frame duration & frame
rate at the chosen `fs`, total pixel/line counts, and the min/max/mean of
both the **intra-line pixel spacing** (dwell time per pixel) and the
**line-turnaround gap** (flyback/settle time) — the two numbers you need to
set a TCSPC's binning window and dead-time tolerance correctly.

---

## 8. How these clocks are actually used with a TCSPC

This section connects the generated files to a real acquisition, and
separates two things that are easy to conflate:

- **Spatial correction** (already solved, upstream of everything here): the
  geometric transform that turned a would-be-distorted square into the
  trapezoidal `lines.csv` voltages, so the *projected* pattern is a uniform
  grid on the tilted wall.
- **Temporal synchronization** (what this whole project builds): making
  sure a photon detected at some macro-time can be attributed to the
  correct one of the 256 scan points, regardless of which of the two scan
  shapes produced it.

### 8.1 Hardware wiring

The three clock columns (`pixel_clk`, `line_clk`, `frame_clk`) are digital
TTL/LVTTL signals generated **on the same sample clock and buffer** that
also drives the galvo's analog `x`/`y` voltages — e.g. a DAQ card's
hardware-correlated AO+DO tasks sharing one `ao/SampleClock`, or an FPGA
BRAM table stepped by one counter that drives the DAC *and* the three
digital outputs in the same clock cycle. This is why the whole pipeline is
organized around a single sample index `n`: it is the one thing every
output — analog voltage and all three digital clocks — is generated
*from*, so there is no separate oscillator or software timer that could
drift relative to the galvo.

Most TCSPC hardware (PicoQuant HydraHarp/MultiHarp, Swabian Time Tagger,
etc.) has dedicated external **marker/sync inputs** for exactly this
purpose. Typical wiring:

```
DAC/FPGA sample counter (n, at fs Hz)
   │
   ├── analog X,Y  ──────────────► galvo driver
   ├── digital pixel_clk ────────► TCSPC "marker 1" / pixel-sync input
   ├── digital line_clk  ────────► TCSPC "marker 2" / line-sync input
   └── digital frame_clk ────────► TCSPC "marker 3" / frame-sync input
```

The TCSPC timestamps every detected photon's arrival (its **macro-time**)
in the *same* clock domain as these markers (either by sharing the sample
clock directly, or by timestamping the marker edges themselves alongside
the photon events in its own high-resolution timebase). That shared
timebase is what makes the next step possible.

### 8.2 Turning marker events into a pixel-resolved image

During acquisition, the TCSPC produces a stream of two kinds of timestamped
events on one timeline: **photon arrivals** and **marker edges**
(pixel/line/frame). Reconstruction is a streaming state machine:

1. On a **frame_clk** edge → reset the pixel/line counters, start a new frame buffer.
2. On a **line_clk** edge → advance the line counter, reset the pixel-within-line counter.
3. On a **pixel_clk** edge → advance the pixel counter; every photon macro-time from *this* edge up to (but not including) the *next* `pixel_clk` edge belongs to the current `(line_index, pixel_index)` — this window is exactly the `pixel_valid_gate` column.
4. Every photon between two `pixel_clk` edges gets binned into that pixel's histogram/counter. Its position on the wall is looked up from the corresponding row of `sync_pixel_edges.csv` (or, for NLOS specifically, from whatever wall-space calibration maps that scan point after the geometric correction — separate from this timing pipeline).

Because `pixel_index` and `line_index` are already computed and stored in
`sync_clocks_full.csv`, a host-side or FPGA-side accumulator does not need
to re-derive scan geometry at acquisition time at all — it only needs to
count marker edges and bin photons between them, using the pre-computed
table as its map from "which edge number" to "which physical scan point."

### 8.3 Why frame_clk is anchored to the first/last *pixel*, not the buffer edges

`frame_clk` fires at the sample of the first `pixel_clk` and the sample of
the last `pixel_clk` (§4b/§7), not at `n = 0` / `n = N_total − 1`. This
matters for TCSPC integration specifically: any lead-in ramp or trailing
settle time in the raw DAC buffer (before the galvo has reached the first
real scan point, or after it has left the last one) would otherwise be
counted as "inside the frame" by a naive frame-boundary detector, wasting
acquisition time on photons that can't be assigned to any pixel. Anchoring
to real pixel events means the TCSPC's frame window exactly brackets valid
data.

### 8.4 What differs between the square and trapezoidal (root_pair) cases

Nothing about *this* protocol changes between the two patterns — the
TCSPC-side state machine in §8.2 is identical either way, which is the
point of building clocks as sample-indexed lookup tables rather than
fixed-frequency counters (see §5). What differs is only the **content** of
the lookup table:

- **Square pattern:** `pixel_clk` pulses are evenly spaced in time
  throughout the frame — a fixed pixel-clock frequency (see
  `clk_out_square/summary.txt`).
- **Trapezoidal (root_pair) pattern:** `pixel_clk` frequency **increases**
  line by line (33.3 Hz on line 0 up to 66.7 Hz on line 15 in this
  dataset's ground truth — see `clk_out_root_groundtruth/line_shrinking_proof.png`),
  because narrower lines near the compressed edge of the trapezoid are
  scanned faster to keep the same constant galvo velocity. A TCSPC
  acquisition doesn't need to know this in advance — it just follows
  whatever `pixel_clk` edges actually arrive — but it does mean **dwell
  time per pixel is not constant across the frame** for this pattern, which
  matters if you're setting a fixed integration/gate time per pixel
  downstream rather than using the marker-defined `pixel_valid_gate`
  window directly.

### 8.5 Timing-budget checklist for TCSPC integration

Using this project's numbers as a concrete example (fs = 1000 Hz):

| Quantity | Value | Where it comes from |
|---|---|---|
| Shortest pixel dwell (root_pair) | 15 ms → 225/15 samples at higher fs | `sync_out_root/summary.txt`, "Intra-line pixel spacing" |
| Longest flyback/settle gap | up to 254 ms | `sync_out_root/summary.txt`, "Line-turnaround gap" |
| Required marker jitter | ≤ 1% of shortest dwell → ≤ 150 µs at fs = 1000 Hz | rule of thumb from the timing-design analysis |
| Frame rate | ~0.12 Hz at fs = 1000 Hz (8.35 s/frame) | `summary.txt`, scales linearly with `fs` |

If `fs` is increased (faster galvo), every duration above scales down
linearly (`t = n / fs`), and the marker-jitter budget tightens
proportionally — at that point, generate the clocks in FPGA/hardware-timed
DAQ rather than software, per the jitter comparison in the original timing
analysis earlier in this project.

---

## 9. File map

| Path | What it is |
|---|---|
| `generate_sync_clocks.py` | Pipeline A: exact clocks from a real `meas_pts.csv` |
| `equal_spacing_from_trajectory.py` | Pipeline B step 1: pixel positions from trajectory geometry only |
| `clocks_from_equal_spacing.py` | Pipeline B step 2: clock timing (direct or velocity-synthetic) from those positions |
| `plot_scan_patterns.py` | Early utility: overlay measurement points on the raw trajectory (matplotlib, light theme) |
| `sync_out_root/`, `sync_out_square/` | Pipeline A outputs (ground truth) |
| `eq_out_root/`, `eq_out_square/` | Pipeline B step 1 outputs (pixel positions) |
| `clk_out_root/`, `clk_out_square/` | Pipeline B step 2 outputs (clock timeline + trajectory-with-events) |
| `clk_out_root_groundtruth/` | Pipeline A outputs re-rendered in the dark trajectory-with-events style, plus `line_shrinking_proof.png` (explicit burst-duration/frequency-per-line charts) |

Each `sync_out_*` / `clk_out_*` directory contains:
- `sync_clocks_full.csv` — the master per-sample table (`n, t_s, x, y, pixel_clk, pixel_index, line_clk, line_index, frame_clk, ...`)
- `sync_clocks_timeline.png` / `sync_clocks_full.png` — X/Y voltage + all three clocks over the full frame
- `sync_clocks_trajectory.png` — dark-themed trajectory with pixel-clock (red) and line-clock (cyan ring) events overlaid
- `summary.txt` — dwell time, flyback gap, and frame-duration statistics

---

## 10. Quick reproduction

```bash
# Ground truth (exact) -- use whenever meas_pts.csv exists
python generate_sync_clocks.py --lines root_lines.csv --meas root_meas_pts.csv --fs 1000 --out sync_out_root

# Trajectory-only reconstruction -- use when it doesn't
python equal_spacing_from_trajectory.py --lines root_lines.csv --n-lines 16 --points-per-line 16 --out eq_out_root
python clocks_from_equal_spacing.py --lines root_lines.csv --pixels eq_out_root/equal_spaced_pixels.csv --fs 1000 --out clk_out_root --title root_pair
```

`--fs` can be changed freely in either pipeline — every downstream
computation is in sample-index space, so raising the galvo sample rate only
rescales `t_s = n / fs`; none of the pixel/line/frame *logic* changes.
