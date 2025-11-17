#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
from collections import defaultdict

# ───────────────────────────────────────────────
# HANDLE COMMAND-LINE ARGUMENT
# ───────────────────────────────────────────────
if len(sys.argv) != 2:
    print("Usage: python TFBS_TSS.py <input_csv>")
    sys.exit(1)

input_file = sys.argv[1]

if not os.path.exists(input_file):
    print(f"❌ Input file not found: {input_file}")
    sys.exit(1)

print(f"📂 Using input file: {input_file}")

# ───────────────────────────────────────────────
# STEP 1: Load and clean data
# ───────────────────────────────────────────────
df = pd.read_csv(input_file, sep="\t")

# Clean column names
df.columns = df.columns.str.strip().str.replace('\ufeff', '')

# Ensure required columns exist
required_columns = ['TFBS', 'TFBS_phase', 'TFBS_strand']
missing = [col for col in required_columns if col not in df.columns]
if missing:
    print(f"❌ Missing columns: {missing}")
    sys.exit(1)

# Clean relevant columns
df['TFBS'] = df['TFBS'].str.strip()
df['TFBS_strand'] = df['TFBS_strand'].str.strip()
df['TFBS_phase'] = df['TFBS_phase'].str.strip()

# Safely convert TFBS_phase to integer
df['phase_int'] = pd.to_numeric(df['TFBS_phase'].str.extract(r'(\d+)')[0], errors='coerce')
df = df.dropna(subset=['phase_int'])
df['phase_int'] = df['phase_int'].astype(int)

# ───────────────────────────────────────────────
# PART 1: Plot total counts per TF (stripped from TFBS)
# ───────────────────────────────────────────────
plt.rcParams.update({
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10
})
BAR_WIDTH = 0.5

df['TF_Name'] = df['TFBS'].str.replace('_M1', '', regex=False)
tf_counts = df['TF_Name'].value_counts().reset_index()
tf_counts.columns = ['TF_Name', 'Site_Count']
tf_counts = tf_counts.sort_values(by='Site_Count', ascending=False)

plt.figure(figsize=(12, 6))
sns.barplot(data=tf_counts, x='TF_Name', y='Site_Count', palette='Reds_r', width=BAR_WIDTH)
plt.title('Number of Sites per TF (TF_M1)')
plt.xlabel('TF Name')
plt.ylabel('Number of Sites')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('TF_M1_barplot_TSS.png')
plt.close()
print("✅ TF site count bar plot saved as: TF_M1_barplot_TSS.png")

# ───────────────────────────────────────────────
# PART 2: TFBS_phase distribution for each TFBS
# ───────────────────────────────────────────────
plt.rcParams.update({
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12
})

tfbs_counts = df['TFBS'].value_counts()
high_tfbs = tfbs_counts.index.tolist()

output_dir = "TFBS_phase_TSS_barplots"
os.makedirs(output_dir, exist_ok=True)

for tfbs in high_tfbs:
    subset = df[df['TFBS'] == tfbs]
    for strand in ['oppo', 'same']:
        strand_subset = subset[subset['TFBS_strand'] == strand]
        if strand_subset.empty:
            continue

        phase_counts = strand_subset['phase_int'].value_counts().sort_index()

        plt.figure(figsize=(8, 5))
        sns.barplot(x=phase_counts.index, y=phase_counts.values, palette='Blues', width=BAR_WIDTH)
        plt.title(f"{tfbs} - Strand: {strand}")
        plt.xlabel('TFBS Phase to TSS')
        plt.ylabel('Count')
        plt.xticks(rotation=0)
        plt.tight_layout()

        filename = f"{tfbs}_strand-{strand}_TFBS_phase.png".replace("/", "_")
        plt.savefig(os.path.join(output_dir, filename))
        plt.close()

print(f"✅ TFBS_phase bar plots saved in folder: {output_dir}/")

# ───────────────────────────────────────────────
# STEP 2: Count TFBS_phase occurrences by TFBS and strand
# ───────────────────────────────────────────────
all_phases = [f"phase{i}" for i in range(10)]
tfbs_counts_dict = defaultdict(lambda: defaultdict(int))

for _, row in df.iterrows():
    key = (row['TFBS'], row['TFBS_strand'])
    tfbs_counts_dict[key][f"phase{row['phase_int']}"] += 1

# ───────────────────────────────────────────────
# STEP 3: Write counts to individual .out files
# ───────────────────────────────────────────────
output_dir_counts = "TFBS_phase_counts_TSS"
os.makedirs(output_dir_counts, exist_ok=True)

for (tfbs, strand), phase_dict in tfbs_counts_dict.items():
    output_file = f"{tfbs}_{strand}.out".replace("/", "_")
    output_path = os.path.join(output_dir_counts, output_file)

    with open(output_path, "w") as f:
        for phase in all_phases:
            count = phase_dict.get(phase, 0)
            f.write(f"{phase} {count}\n")

print(f"✅ Output files written to: {output_dir_counts}/")

# ───────────────────────────────────────────────
# PART 3: Strand-specific total counts per TFBS
# ───────────────────────────────────────────────
BAR_WIDTH = 0.5
strand_counts = (
    df.groupby(['TFBS', 'TFBS_strand'])
    .size()
    .unstack(fill_value=0)
    .reset_index()
)

for col in ['same', 'oppo']:
    if col not in strand_counts.columns:
        strand_counts[col] = 0

strand_counts['strand_ratio'] = strand_counts[['same', 'oppo']].max(axis=1) / strand_counts[['same', 'oppo']].min(axis=1).replace(0, 1e-6)
strand_counts = strand_counts.sort_values(by='strand_ratio', ascending=False)

plot_df = strand_counts.melt(id_vars='TFBS', value_vars=['same', 'oppo'], var_name='Strand', value_name='Count')
plot_df['TFBS'] = pd.Categorical(plot_df['TFBS'], categories=strand_counts['TFBS'], ordered=True)

plt.figure(figsize=(14, 6))
sns.barplot(data=plot_df, x='TFBS', y='Count', hue='Strand', palette='Set2', width=BAR_WIDTH)
plt.title('TFBS Strand Bias (Grouped Counts, Sorted by Imbalance)')
plt.xlabel('TFBS')
plt.ylabel('Site Count')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('TFBS_strand_bias_ratio_barplot.png')
plt.close()
print("✅ Strand bias grouped bar plot saved as: TFBS_strand_bias_ratio_barplot.png")

# ───────────────────────────────────────────────
# PART 4: Max phase skew (5-phase vs rest), by strand
# ───────────────────────────────────────────────
phase_skew_results = []

for tfbs in high_tfbs:
    for strand in ['same', 'oppo']:
        subset = df[(df['TFBS'] == tfbs) & (df['TFBS_strand'] == strand)]
        if subset.empty:
            continue

        phase_counts = subset['phase_int'].value_counts().reindex(range(10), fill_value=0).to_list()

        max_ratio = 0
        best_window = ()

        for start in range(10):
            window_phases = [(start + i) % 10 for i in range(5)]
            window_sum = sum(phase_counts[p] for p in window_phases)
            rest_sum = sum(phase_counts) - window_sum
            ratio = window_sum / rest_sum if rest_sum > 0 else float('inf')
            if ratio > max_ratio:
                max_ratio = ratio
                best_window = window_phases

        phase_skew_results.append({
            'TFBS': tfbs,
            'Strand': strand,
            'Max_Phase_Ratio': max_ratio,
            'Best_Window': '-'.join(str(p) for p in best_window)
        })

skew_df = pd.DataFrame(phase_skew_results)

# Plot same strand
same_df = skew_df[skew_df['Strand'] == 'same'].sort_values(by='Max_Phase_Ratio', ascending=False)
same_df['TFBS'] = pd.Categorical(same_df['TFBS'], categories=same_df['TFBS'], ordered=True)

plt.figure(figsize=(12, 5))
sns.barplot(data=same_df, x='TFBS', y='Max_Phase_Ratio', palette='Blues', width=BAR_WIDTH)
plt.title('Max 5-Phase Bias (TFBS in same strand)')
plt.xlabel('TFBS')
plt.ylabel('Max 5-Phase / Rest Ratio')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('TFBS_phase_skew_barplot_same_TSS.png')
plt.close()

# Plot opposite strand
oppo_df = skew_df[skew_df['Strand'] == 'oppo'].sort_values(by='Max_Phase_Ratio', ascending=False)
oppo_df['TFBS'] = pd.Categorical(oppo_df['TFBS'], categories=oppo_df['TFBS'], ordered=True)

plt.figure(figsize=(12, 5))
sns.barplot(data=oppo_df, x='TFBS', y='Max_Phase_Ratio', palette='Reds', width=BAR_WIDTH)
plt.title('Max 5-Phase Bias (TFBS in opposite strand)')
plt.xlabel('TFBS')
plt.ylabel('Max 5-Phase / Rest Ratio')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('TFBS_phase_skew_barplot_oppo_TSS.png')
plt.close()

print("✅ Saved:")
print("   TFBS_phase_skew_barplot_same_TSS.png")
print("   TFBS_phase_skew_barplot_oppo_TSS.png")

# ───────────────────────────────────────────────
# FINAL STEP: Save ratio values
# ───────────────────────────────────────────────
strand_counts[['TFBS', 'same', 'oppo', 'strand_ratio']].to_csv(
    'TFBS_strand_bias_ratios.tsv', sep='\t', index=False
)
skew_df[['TFBS', 'Strand', 'Max_Phase_Ratio', 'Best_Window']].to_csv(
    'TFBS_phase_skew_ratios_TSS.tsv', sep='\t', index=False
)

print("📄 Saved ratio values:")
print("   → TFBS_strand_bias_ratios.tsv")
print("   → TFBS_phase_skew_ratios_TSS.tsv")
