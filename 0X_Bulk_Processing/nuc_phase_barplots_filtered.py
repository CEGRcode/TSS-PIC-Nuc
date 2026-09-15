#!/usr/bin/env python3
"""
nuc_phase_barplots_filtered.py
================================

Same set of plots as nuc_phase_barplots.py, EXCEPT it only includes
rows where the Nuc number in column $5 satisfies:  40 <= $5 <= 160

All downstream plots (phase score, TSS signal, and entry counts --
each broken down by Nuc type / coding-gene info / HMM annotation, plus
an ungrouped "all data" plot) are generated only from this filtered
subset of rows.

Scans a folder for .bed files and, for EACH category within each of the
three groupings below, generates a SEPARATE bar plot with 10 bars (one
per phase), shown in the order: phase9, phase0, phase1, ..., phase8.

Two metrics are plotted (each gets its own full set of plots):
  A. Nuc phase score  ($7)
  B. TSS signal       ($18)

In addition, a third set of plots shows the raw NUMBER OF ENTRIES
(count) falling into each phase bin, per category (metric-independent).

Groupings (each category gets its own plot), plus one ungrouped
"overall" plot per metric:
  1. Nuc type (column 8):
        sameYR_lowWS, lowYR_sameWS, YRWS, antiYR_antiWS,
        antiYR_sameWS, sameYR_antiWS   (labeled Type 1 - Type 6)
  2. Coding-gene information (column 27): codingTSS, other
  3. HMM annotation (column 27): Enhancer, Promoter

Column reference (1-based, as given by the user / awk-style $N):
  $5  -> Nuc number, used ONLY as a filter: keep 40 <= $5 <= 160
         (0-based index 4)
  $7  -> Nuc phase score      (0-based index 6)
  $8  -> Nuc type              (0-based index 7)
  $18 -> TSS signal            (0-based index 17)
  $21 -> phase information, contains the substring "phase0".."phase9"
         (0-based index 20)
  $27 -> coding-gene info ("codingTSS"/"other") AND
         HMM annotation ("Enhancer"/"Promoter") -- both live in the
         same field, distinguished by which keyword is present
         (0-based index 26)

For a given category (e.g. Nuc type == "YRWS"), each bar is the MEAN
of the metric (phase score or TSS signal) for all (filtered) rows in
that category that fall in that phase bin, with an error bar (SEM by
default, or SD with --error sd) computed from those same raw values.
The x-axis order is fixed as:
    phase9, phase0, phase1, phase2, phase3, phase4, phase5, phase6,
    phase7, phase8

Usage
-----
    python nuc_phase_barplots_filtered.py /path/to/folder [-o OUTPUT_DIR] [--error sem|sd]

Requirements: pandas, numpy, matplotlib  (all pure-Python, no seaborn
needed)
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Column indices (0-based) matching the user's 1-based $N notation
# ----------------------------------------------------------------------
COL_NUC_NUMBER = 4     # $5  (used only for filtering: 40 <= value <= 160)
COL_PHASE_SCORE = 6    # $7
COL_NUC_TYPE = 7       # $8
COL_TSS_SIGNAL = 17    # $18
COL_PHASE_INFO = 20    # $21
COL_ANNOT = 26         # $27

# Inclusive Nuc-number filter range
NUC_NUMBER_MIN = 40
NUC_NUMBER_MAX = 160

# Minimum number of columns a line must have to be usable
MIN_COLS = COL_ANNOT + 1

NUC_TYPES = [
    "sameYR_lowWS",
    "lowYR_sameWS",
    "YRWS",
    "antiYR_antiWS",
    "antiYR_sameWS",
    "sameYR_antiWS",
]

# Display labels: "Type N: <nuc_type>" in the order given by the user
NUC_TYPE_LABELS = {
    nuc_type: f"Type {i}: {nuc_type}" for i, nuc_type in enumerate(NUC_TYPES, start=1)
}

PHASE_RE = re.compile(r"phase\s*([0-9])", re.IGNORECASE)

# Natural order (used internally / for CSV export)
PHASE_ORDER = [f"phase{i}" for i in range(10)]

# Display order requested: phase9, phase0, phase1, ..., phase8
PHASE_DISPLAY_ORDER = ["phase9"] + [f"phase{i}" for i in range(0, 9)]


def load_bed_folder(folder):
    """Read every .bed file in `folder`, extract the needed columns,
    and return one concatenated DataFrame with clean columns:
    ['phase_score', 'nuc_type', 'phase', 'annot']
    """
    bed_files = sorted(glob.glob(os.path.join(folder, "*.bed")))
    if not bed_files:
        sys.exit(f"No .bed files found in: {folder}")

    frames = []
    for fp in bed_files:
        try:
            df = pd.read_csv(
                fp,
                sep="\t",
                header=None,
                comment="#",
                engine="python",
                on_bad_lines="skip",
                dtype=str,
            )
        except Exception as e:
            print(f"  [WARN] could not read {fp}: {e}", file=sys.stderr)
            continue

        if df.shape[1] < MIN_COLS:
            print(
                f"  [WARN] {fp} has only {df.shape[1]} columns "
                f"(need >= {MIN_COLS}); skipped",
                file=sys.stderr,
            )
            continue

        sub = pd.DataFrame(
            {
                "nuc_number": pd.to_numeric(
                    df.iloc[:, COL_NUC_NUMBER], errors="coerce"
                ),
                "phase_score": pd.to_numeric(
                    df.iloc[:, COL_PHASE_SCORE], errors="coerce"
                ),
                "tss_signal": pd.to_numeric(
                    df.iloc[:, COL_TSS_SIGNAL], errors="coerce"
                ),
                "nuc_type": df.iloc[:, COL_NUC_TYPE].astype(str).str.strip(),
                "phase_raw": df.iloc[:, COL_PHASE_INFO].astype(str),
                "annot": df.iloc[:, COL_ANNOT].astype(str),
            }
        )
        sub["source_file"] = os.path.basename(fp)
        frames.append(sub)

    if not frames:
        sys.exit("No usable data could be parsed from any .bed file.")

    data = pd.concat(frames, ignore_index=True)

    # Extract phase0..phase9 from the phase_raw field
    data["phase"] = data["phase_raw"].apply(_extract_phase)

    # Drop rows with missing phase info (needed for both metrics).
    # Missing individual metric values (phase_score / tss_signal) are
    # left as NaN and simply skipped by the per-metric aggregations.
    data = data.dropna(subset=["phase"])

    # --------------------------------------------------------------
    # Apply the Nuc-number filter: keep only 40 <= $5 <= 160
    # --------------------------------------------------------------
    n_before = len(data)
    data = data[
        data["nuc_number"].notna()
        & (data["nuc_number"] >= NUC_NUMBER_MIN)
        & (data["nuc_number"] <= NUC_NUMBER_MAX)
    ]
    n_after = len(data)
    print(
        f"Nuc-number filter ({NUC_NUMBER_MIN} <= $5 <= {NUC_NUMBER_MAX}): "
        f"kept {n_after} / {n_before} rows"
    )

    return data


def _extract_phase(text):
    m = PHASE_RE.search(text)
    if m:
        return f"phase{m.group(1)}"
    return None


def phase_counts_for_category(data, group_col, category):
    """For a single category, count how many rows fall in each phase
    bin (regardless of whether phase_score/tss_signal happen to be
    NaN for a given row -- this is a raw entry count per phase).

    Returns a Series indexed by phase (in PHASE_DISPLAY_ORDER).
    """
    sub = data[data[group_col] == category]
    counts = sub.groupby("phase").size()
    counts = counts.reindex(PHASE_DISPLAY_ORDER).fillna(0).astype(int)
    return counts


def make_count_barplot(counts, title, out_path):
    """Bar plot for ONE category: x-axis = the 10 phases (in
    PHASE_DISPLAY_ORDER), y-axis = number of entries (count) in that
    phase bin. No error bars (a count has no variance to show); the
    count itself is annotated above each bar.
    """
    labels = [p.replace("phase", "P") for p in counts.index]
    values = counts.values

    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.bar(
        x,
        values,
        color="#ffb703",
        edgecolor="black",
        alpha=0.85,
        zorder=2,
    )

    ax.set_title(title, fontsize=13)
    ax.set_xlabel("Phase", fontsize=11)
    ax.set_ylabel("Number of entries", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)

    ymax = values.max() if len(values) else 0
    ax.set_ylim(0, ymax * 1.15 if ymax > 0 else 1)
    for i, v in enumerate(values):
        ax.text(
            i, v, f"{v}", ha="center", va="bottom", fontsize=9, color="black",
        )

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"  saved: {out_path}")


def phase_stats_overall(data, value_col, error="sem"):
    """Same as phase_stats_for_category but computed over ALL rows,
    with no grouping/category filter applied.
    """
    grouped = data.groupby("phase")[value_col]
    means = grouped.mean()
    counts = grouped.count()
    if error == "sd":
        errs = grouped.std(ddof=1)
    else:
        errs = grouped.std(ddof=1) / np.sqrt(counts)

    table = pd.DataFrame({"mean": means, "error": errs, "n": counts})
    table = table.reindex(PHASE_DISPLAY_ORDER)
    return table


def phase_counts_overall(data):
    """Same as phase_counts_for_category but computed over ALL rows."""
    counts = data.groupby("phase").size()
    counts = counts.reindex(PHASE_DISPLAY_ORDER).fillna(0).astype(int)
    return counts


def phase_stats_for_category(data, group_col, category, value_col, error="sem"):
    """For a single category (e.g. Nuc type == 'YRWS'), compute, for
    each phase bin (phase0..phase9), the mean of `value_col` and an
    error value (SEM or SD) across all raw rows of that category/phase.

    Returns a DataFrame indexed by phase (in PHASE_DISPLAY_ORDER) with
    columns: mean, error, n
    """
    sub = data[data[group_col] == category]

    grouped = sub.groupby("phase")[value_col]
    means = grouped.mean()
    counts = grouped.count()
    if error == "sd":
        errs = grouped.std(ddof=1)
    else:
        errs = grouped.std(ddof=1) / np.sqrt(counts)

    table = pd.DataFrame({"mean": means, "error": errs, "n": counts})
    table = table.reindex(PHASE_DISPLAY_ORDER)
    return table


def make_category_barplot(stats_table, title, ylabel, out_path, error_label):
    """Bar plot for ONE category: x-axis = the 10 phases (in
    PHASE_DISPLAY_ORDER), y-axis = mean phase score, with error bars
    (SEM or SD, per `error_label`) and n= annotations per bar.
    """
    labels = [p.replace("phase", "P") for p in stats_table.index]
    means = stats_table["mean"].values
    errs = stats_table["error"].fillna(0).values
    ns = stats_table["n"].fillna(0).astype(int).values

    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.bar(
        x,
        means,
        yerr=errs,
        capsize=4,
        color="#8ecae6",
        edgecolor="black",
        alpha=0.85,
        zorder=2,
    )

    ax.set_title(title, fontsize=13)
    ax.set_xlabel("Phase", fontsize=11)
    ax.set_ylabel(f"{ylabel}\n(bar = mean \u00b1 {error_label})", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)

    valid_means = means[~np.isnan(means)]
    if len(valid_means) == 0:
        ax.text(
            0.5, 0.5, "No data", ha="center", va="center",
            transform=ax.transAxes, fontsize=12, color="gray",
        )
    else:
        ymax = np.nanmax(means + np.nan_to_num(errs))
        ymin = min(0, np.nanmin(means - np.nan_to_num(errs)))
        span = (ymax - ymin) if ymax > ymin else 1
        ax.set_ylim(ymin - 0.05 * span, ymax + 0.15 * span)
        for i, n in enumerate(ns):
            ax.text(
                i,
                ax.get_ylim()[1],
                f"n={n}",
                ha="center",
                va="top",
                fontsize=8,
                color="gray",
            )

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"  saved: {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder", help="Folder containing .bed files")
    ap.add_argument(
        "-o",
        "--outdir",
        default=None,
        help="Output directory for plots/tables (default: <folder>/phase_barplots_filtered)",
    )
    ap.add_argument(
        "--error",
        choices=["sem", "sd"],
        default="sem",
        help="Error bar type: 'sem' (default) or 'sd'",
    )
    args = ap.parse_args()

    folder = args.folder
    outdir = args.outdir or os.path.join(folder, "phase_barplots_filtered")
    os.makedirs(outdir, exist_ok=True)

    print(f"Loading .bed files from: {folder}")
    data = load_bed_folder(folder)
    print(f"Total usable rows parsed: {len(data)}")

    error_label = "SD" if args.error == "sd" else "SEM"

    def run_group(value_col, ylabel, metric_outdir, group_col, categories,
                  subfolder, title_prefix, label_map=None):
        group_outdir = os.path.join(metric_outdir, subfolder)
        os.makedirs(group_outdir, exist_ok=True)
        for cat in categories:
            table = phase_stats_for_category(
                data, group_col, cat, value_col, error=args.error
            )
            display_name = (label_map or {}).get(cat, cat)
            safe_name = re.sub(r"[^A-Za-z0-9_.-]", "_", str(display_name))
            table.to_csv(os.path.join(group_outdir, f"{safe_name}_phase_stats.csv"))
            make_category_barplot(
                table,
                f"{title_prefix}: {display_name}" if label_map is None else display_name,
                ylabel,
                os.path.join(group_outdir, f"barplot_{safe_name}.png"),
                error_label,
            )

    # Prepare the two extra category columns once (shared by both metrics)
    data["gene_info"] = np.select(
        [
            data["annot"].str.contains("codingTSS", na=False),
            data["annot"].str.contains("other", na=False),
        ],
        ["codingTSS", "other"],
        default=None,
    )
    data["hmm_info"] = np.select(
        [
            data["annot"].str.contains("Enhancer", na=False),
            data["annot"].str.contains("Promoter", na=False),
        ],
        ["Enhancer", "Promoter"],
        default=None,
    )

    metrics = [
        ("phase_score", "Mean phase score", "phase_score"),
        ("tss_signal", "Mean TSS signal", "tss_signal"),
    ]

    def run_count_group(group_col, categories, subfolder, title_prefix, label_map=None):
        group_outdir = os.path.join(outdir, "count", subfolder)
        os.makedirs(group_outdir, exist_ok=True)
        for cat in categories:
            counts = phase_counts_for_category(data, group_col, cat)
            display_name = (label_map or {}).get(cat, cat)
            safe_name = re.sub(r"[^A-Za-z0-9_.-]", "_", str(display_name))
            counts.rename("count").to_csv(
                os.path.join(group_outdir, f"{safe_name}_phase_counts.csv")
            )
            make_count_barplot(
                counts,
                f"{title_prefix}: {display_name}" if label_map is None else display_name,
                os.path.join(group_outdir, f"barplot_{safe_name}.png"),
            )

    for value_col, ylabel, metric_folder in metrics:
        metric_outdir = os.path.join(outdir, metric_folder)
        os.makedirs(metric_outdir, exist_ok=True)
        print(f"\n--- Metric: {value_col} ---")

        # 0) One overall bar plot (no grouping at all)
        overall_dir = os.path.join(metric_outdir, "overall")
        os.makedirs(overall_dir, exist_ok=True)
        overall_table = phase_stats_overall(data, value_col, error=args.error)
        overall_table.to_csv(os.path.join(overall_dir, "all_phase_stats.csv"))
        make_category_barplot(
            overall_table,
            "All data (no grouping)",
            ylabel,
            os.path.join(overall_dir, "barplot_all.png"),
            error_label,
        )

        # 1) One bar plot per Nuc type (labeled Type 1, Type 2, ...)
        run_group(
            value_col, ylabel, metric_outdir,
            "nuc_type", NUC_TYPES, "nuc_type", "Nuc type",
            label_map=NUC_TYPE_LABELS,
        )

        # 2) One bar plot per coding-gene category (codingTSS, other)
        run_group(
            value_col, ylabel, metric_outdir,
            "gene_info", ["codingTSS", "other"], "coding_gene", "Coding gene info",
        )

        # 3) One bar plot per HMM annotation category (Enhancer, Promoter)
        run_group(
            value_col, ylabel, metric_outdir,
            "hmm_info", ["Enhancer", "Promoter"], "hmm_annotation", "HMM annotation",
        )

    # ------------------------------------------------------------
    # 4) Number of entries per phase bin, one plot per category
    #    (this is metric-independent -- just a raw count)
    # ------------------------------------------------------------
    print("\n--- Metric: count (number of entries per phase) ---")

    # Overall count plot (no grouping)
    overall_count_dir = os.path.join(outdir, "count", "overall")
    os.makedirs(overall_count_dir, exist_ok=True)
    overall_counts = phase_counts_overall(data)
    overall_counts.rename("count").to_csv(
        os.path.join(overall_count_dir, "all_phase_counts.csv")
    )
    make_count_barplot(
        overall_counts,
        "All data (no grouping)",
        os.path.join(overall_count_dir, "barplot_all.png"),
    )

    run_count_group("nuc_type", NUC_TYPES, "nuc_type", "Nuc type", label_map=NUC_TYPE_LABELS)
    run_count_group("gene_info", ["codingTSS", "other"], "coding_gene", "Coding gene info")
    run_count_group("hmm_info", ["Enhancer", "Promoter"], "hmm_annotation", "HMM annotation")

    print(f"\nDone. All outputs written to: {outdir}")


if __name__ == "__main__":
    main()
