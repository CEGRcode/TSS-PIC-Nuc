import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


# ───────────────────────────────────────────────
# STEP 1: Load and clean data
# ───────────────────────────────────────────────
input_file = 'TF_M1_adjNuc.csv'
df = pd.read_csv(input_file)

# Clean column names
df.columns = df.columns.str.strip().str.replace('\ufeff', '')

# Ensure required columns exist
required_columns = ['TFBS', 'TFBS_phase', 'TFBS_strand']
missing = [col for col in required_columns if col not in df.columns]
if missing:
    print(f"❌ Missing columns: {missing}")
    raise SystemExit

# Clean strand values
df['TFBS_strand'] = df['TFBS_strand'].str.strip()

# ───────────────────────────────────────────────
# PART 1: Plot total counts per TF (stripped from TFBS)
# ───────────────────────────────────────────────

# ───────────────────────────────────────────────
# GLOBAL SETTINGS: Fonts and Bar Width
# ───────────────────────────────────────────────
plt.rcParams.update({
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10
})
BAR_WIDTH = 0.5  # Narrower bars


df['TF_Name'] = df['TFBS'].str.replace('_M1', '', regex=False)
tf_counts = df['TF_Name'].value_counts().reset_index()
tf_counts.columns = ['TF_Name', 'Site_Count']
tf_counts = tf_counts.sort_values(by='Site_Count', ascending=False)

# Plot
plt.figure(figsize=(12, 6))
sns.barplot(data=tf_counts, x='TF_Name', y='Site_Count', palette='Reds_r', width=BAR_WIDTH)
plt.title('Number of Sites per TF (TF_M1)')
plt.xlabel('TF Name')
plt.ylabel('Number of Sites')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('TF_M1_barplot_nuc.png')
plt.close()
print("✅ TF site count bar plot saved as: TF_M1_barplot.png")

# ───────────────────────────────────────────────
# PART 2: TFBS_phase distribution for high-count TFBS
# ───────────────────────────────────────────────

# ───────────────────────────────────────────────
# GLOBAL SETTINGS: Fonts and Bar Width
# ───────────────────────────────────────────────
plt.rcParams.update({
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12
})
BAR_WIDTH = 0.5  # Narrower bars

# Filter TFBS with >500 counts
tfbs_counts = df['TFBS'].value_counts()
high_tfbs = tfbs_counts[tfbs_counts > 0].index.tolist()

# Create output directory
output_dir = "TFBS_phase_nuc_barplots"
os.makedirs(output_dir, exist_ok=True)

for tfbs in high_tfbs:
    subset = df[df['TFBS'] == tfbs]

    for strand in ['oppo', 'same']:
        strand_subset = subset[subset['TFBS_strand'] == strand]

        if strand_subset.empty:
            continue

        # Convert 'TFBS_phase' from "phaseN" to integer N
        phase_numeric = strand_subset['TFBS_phase'].str.replace("phase", "", regex=False).astype(int)

        # Count and sort
        phase_counts = phase_numeric.value_counts().sort_index()

        # Plot
        plt.figure(figsize=(8, 5))
        sns.barplot(x=phase_counts.index, y=phase_counts.values, palette='Blues_d', width=BAR_WIDTH)
        plt.title(f"{tfbs} - Strand: {strand}")
        plt.xlabel('TFBS Phase to Nuc')
        plt.ylabel('Count')
        plt.xticks(rotation=0)
        plt.tight_layout()

        # Save
        filename = f"{tfbs}_strand-{strand}_TFBS_phase.png".replace("/", "_")
        plt.savefig(os.path.join(output_dir, filename))
        plt.close()

print(f"✅ TFBS_phase bar plots saved in folder: {output_dir}/")

# ───────────────────────────────────────────────
# PART 3: Strand-specific total counts per TFBS
# ───────────────────────────────────────────────
# ───────────────────────────────────────────────
# GLOBAL SETTINGS: Fonts and Bar Width
# ───────────────────────────────────────────────
plt.rcParams.update({
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12
})
BAR_WIDTH = 0.5  # Narrower bars

# Count total per TFBS
total_counts = df['TFBS'].value_counts()
high_tfbs = total_counts[total_counts > 0].index.tolist()

# Filter data
filtered_df = df[df['TFBS'].isin(high_tfbs)]

# Group by TFBS and strand
strand_counts = (
    filtered_df.groupby(['TFBS', 'TFBS_strand'])
    .size()
    .unstack(fill_value=0)
    .reset_index()
)

# Ensure both 'same' and 'oppo' columns exist
if 'same' not in strand_counts.columns:
    strand_counts['same'] = 0
if 'oppo' not in strand_counts.columns:
    strand_counts['oppo'] = 0

# Calculate max/min ratio (strand imbalance)
strand_counts['strand_ratio'] = strand_counts[['same', 'oppo']].max(axis=1) / strand_counts[['same', 'oppo']].min(axis=1).replace(0, 1e-6)

# Sort by strand ratio descending
strand_counts = strand_counts.sort_values(by='strand_ratio', ascending=False)

# Melt for plotting
plot_df = strand_counts.melt(id_vars='TFBS', value_vars=['same', 'oppo'], var_name='Strand', value_name='Count')

# Sort TFBS axis
plot_df['TFBS'] = pd.Categorical(plot_df['TFBS'], categories=strand_counts['TFBS'], ordered=True)

# Plot grouped barplot
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
# Circular window logic included (e.g., 9-0-1-2-3)
# Generates two barplots: one for 'same', one for 'oppo'
# ───────────────────────────────────────────────
plt.rcParams.update({
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12
})
from collections import defaultdict

# Filter TFBS with enough data
high_tfbs = df['TFBS'].value_counts()
high_tfbs = high_tfbs[high_tfbs > 0].index.tolist()

# Convert TFBS_phase to integer
df['phase_int'] = df['TFBS_phase'].str.replace('phase', '', regex=False).astype(int)

phase_skew_results = []

for tfbs in high_tfbs:
    for strand in ['same', 'oppo']:
        subset = df[(df['TFBS'] == tfbs) & (df['TFBS_strand'] == strand)]
        if subset.empty:
            continue

        # Count occurrences in each phase (0–9)
        phase_counts = subset['phase_int'].value_counts().reindex(range(10), fill_value=0).to_list()

        # Circular sliding window
        phase_indices = list(range(10)) + list(range(5))  # simulate circular wrap-around
        max_ratio = 0
        best_window = ()

        for start in range(10):
            window_phases = [phase_indices[start + i] % 10 for i in range(5)]
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

# Convert results to DataFrame
skew_df = pd.DataFrame(phase_skew_results)

# ───────────────────────────────────────────────
# Barplot: same strand
# ───────────────────────────────────────────────
same_df = skew_df[skew_df['Strand'] == 'same'].sort_values(by='Max_Phase_Ratio', ascending=False)
same_df['TFBS'] = pd.Categorical(same_df['TFBS'], categories=same_df['TFBS'], ordered=True)

plt.figure(figsize=(12, 5))
sns.barplot(data=same_df, x='TFBS', y='Max_Phase_Ratio', palette='Blues_d', width=BAR_WIDTH)
plt.title('Max 5-Phase Bias (TFBS in same strand)')
plt.xlabel('TFBS')
plt.ylabel('Max 5-Phase / Rest Ratio')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('TFBS_phase_skew_barplot_same_Nuc.png')
plt.close()

# ───────────────────────────────────────────────
# Barplot: opposite strand
# ───────────────────────────────────────────────
oppo_df = skew_df[skew_df['Strand'] == 'oppo'].sort_values(by='Max_Phase_Ratio', ascending=False)
oppo_df['TFBS'] = pd.Categorical(oppo_df['TFBS'], categories=oppo_df['TFBS'], ordered=True)

plt.figure(figsize=(12, 5))
sns.barplot(data=oppo_df, x='TFBS', y='Max_Phase_Ratio', palette='Reds_d', width=BAR_WIDTH)
plt.title('Max 5-Phase Bias (TFBS in opposite strand)')
plt.xlabel('TFBS')
plt.ylabel('Max 5-Phase / Rest Ratio')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('TFBS_phase_skew_barplot_oppo_Nuc.png')
plt.close()

print("✅ Saved:")
print("   TFBS_phase_skew_barplot_same.png")
print("   TFBS_phase_skew_barplot_oppo.png")


# ───────────────────────────────────────────────
# FINAL STEP: Save ratio values used in Part 5 and 6
# ───────────────────────────────────────────────

# Save Part 5 strand bias ratios
strand_counts[['TFBS', 'same', 'oppo', 'strand_ratio']].to_csv(
    'TFBS_strand_bias_ratios.tsv', sep='\t', index=False
)

# Save Part 6 max 5-phase skew ratios
skew_df[['TFBS', 'Strand', 'Max_Phase_Ratio', 'Best_Window']].to_csv(
    'TFBS_phase_skew_ratios_Nuc.tsv', sep='\t', index=False
)

print("📄 Saved ratio values:")
print("   → TFBS_strand_bias_ratios_Nuc.tsv (Part 5)")
print("   → TFBS_phase_skew_ratios_Nuc.tsv (Part 6)")


