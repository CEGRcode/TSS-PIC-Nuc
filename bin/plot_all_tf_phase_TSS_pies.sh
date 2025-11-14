#!/bin/bash

# Usage check
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <TFBS_phase_skew_ratios.tsv> <TFBS_phase_counts_dir>"
    exit 1
fi

# Command-line arguments
TF_TSS_RATIO="$1"
PHASE_COUNT_DIR="$2"

# Check input file
if [ ! -f "$TF_TSS_RATIO" ]; then
  echo "Error: $TF_TSS_RATIO not found!"
  exit 1
fi

# Check input directory
if [ ! -d "$PHASE_COUNT_DIR" ]; then
  echo "Error: Phase count directory not found: $PHASE_COUNT_DIR"
  exit 1
fi

mkdir -p TF-TSS_pie_charts  # Output directory

# Loop through each line (skip header)
tail -n +2 "$TF_TSS_RATIO" | while IFS=$'\t' read -r TF orientation _ best_window; do

  # Construct phase count file path
  PHASE_COUNT_FILE="${PHASE_COUNT_DIR}/${TF}_${orientation}.out"

  if [ ! -f "$PHASE_COUNT_FILE" ]; then
    echo "Warning: File not found: $PHASE_COUNT_FILE"
    continue
  fi

  # Convert best_window (e.g., "3-4-5-6-7") into Python list: [3, 4, 5, 6, 7]
  BEST_PHASES_LIST=$(echo "$best_window" | awk -F'-' '{printf "["; for (i=1; i<=NF; i++) printf "%s%s", $i, (i<NF?", ":""); print "]"}')

  # Run Python inline
  python3 <<EOF
import matplotlib.pyplot as plt
from itertools import combinations

# Best 5 phases
best_phases = $BEST_PHASES_LIST
best_phases = [int(p) % 10 for p in best_phases]

# Read phase counts
counts = {}
with open("${PHASE_COUNT_FILE}") as f:
    for line in f:
        phase, count = line.strip().split()
        counts[int(phase[-1])] = int(count)

# Phase list and values
all_phases = list(range(10))
values = [counts[p] for p in all_phases]

# Observed best window total
best_sum = sum(counts[p] for p in best_phases)

# Exhaustive permutation test: all 5-phase combinations
comb_sums = [sum(values[i] for i in comb) for comb in combinations(all_phases, 5)]
p_value_perm = sum(1 for s in comb_sums if s >= best_sum) / len(comb_sums)

# Save p-value to file
with open(f"TF-TSS_pie_charts/${TF}_${orientation}_pvalue.out", "w") as out_f:
    out_f.write(f"Best_Window_sum: {best_sum}\\n")
    out_f.write(f"Permutation_p_value: {p_value_perm:.6g}\\n")

# Pie chart with transparent background and black border
colors = ['yellow' if p in best_phases else 'black' for p in range(10)]

fig, ax = plt.subplots(figsize=(6, 6), facecolor='none')  # Transparent background
wedges, _ = ax.pie(values, colors=colors, startangle=90)

# Set thicker edges
for i, wedge in enumerate(wedges):
    wedge.set_edgecolor('black' if colors[i] == 'yellow' else 'white')
    wedge.set_linewidth(3)

# Add black circle as border
circle = plt.Circle((0, 0), 1, color='black', fill=False, linewidth=4)
ax.add_artist(circle)

ax.set(aspect="equal")
plt.tight_layout()
plt.savefig(f"TF-TSS_pie_charts/${TF}_${orientation}_phase_pie.png", dpi=300, transparent=True)
print("✔️ Pie + permutation p-value saved: ${TF}_${orientation}")

EOF

done
