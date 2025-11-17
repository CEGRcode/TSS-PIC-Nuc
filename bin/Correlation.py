#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency
import argparse

# ===============================
# Command-line arguments
# ===============================
parser = argparse.ArgumentParser(description="Enrichment analysis of Nuc_type vs categorical features")
parser.add_argument("-i", "--input", required=True, help="Input CSV file")
args = parser.parse_args()


# ===============================
# Load and filter the data
# ===============================
df = pd.read_csv(args.input, sep='\t')
df.columns = df.columns.str.strip()

# ===============================
# Merge 'NonfixedTATA*' into 'TATA'
# ===============================
tata_groups = [
    'TATAoppo2mis', 'TATAoppo1mis', 'TATAoppo0mis',
    'TATAosame1mis', 'TATAosame2mis', 'TATAosame3mis',
    'NonfixedTATAoppo2mis', 'NonfixedTATAoppo1mis', 'NonfixedTATAoppo0mis',
    'NonfixedTATAosame1mis', 'NonfixedTATAosame2mis', 'NonfixedTATAosame3mis'
]
if 'core-promoter' in df.columns:
    df['core-promoter'] = df['core-promoter'].replace(tata_groups, 'TATA')

# ===============================
# Numerical features and ordering
# ===============================
numerical_features = ['Distance_to_TSS', 'Phase_score', 'Conservation','PropT', 'Expression']

# Set fixed desired Nuc_type order
ordered_nuctypes = ['YRWS', 'sameYR_lowWS', 'lowYR_sameWS', 'sameYR_antiWS', 
                    'antiYR_sameWS', 'antiYR_antiWS', 'lessDNAencode']
ordered_nuctypes = [n for n in ordered_nuctypes if n in df['Nuc_type'].unique()]

# ===============================
# Color palette for Nuc_type
# ===============================
nuc_type_colors = {
    'YRWS': 'lightcoral',            # light red
    'sameYR_lowWS': 'orange',
    'lowYR_sameWS': 'yellow',
    'sameYR_antiWS': 'lightgreen',
    'antiYR_sameWS': 'cyan',
    'antiYR_antiWS': 'lightblue',
    'lessDNAencode': 'grey'
}

# ===============================
# Categorical enrichment heatmaps
# ===============================
df['Orientation'] = df['Orientation'].replace('Secondpair', 'Reference')

# Select categorical features
cat_features = df.select_dtypes(include='object').columns.difference(['Nuc_type', 'Site'])

# Category-specific orders
feature_order_dict = {
    'core-promoter': ['noTATA', 'TATA'],
    'TSS': ['codingTSS', 'other'],
    'HMM': ['Promoter', 'Enhancer'],
    'Orientation': ['Reference', 'Singleton'],
    'core-Nuc': ['NFR', 'NDR'],
    'Annotation1': ['codingTSSCpG', 'othernoCpG']
}

print("\n=== Enrichment analysis of categorical features vs Nuc_type ===\n")

import os

# Open output file for writing
output_file = "enrichment_tables/enrichment_results.out"
with open(output_file, "w") as f_out:

    for feature in cat_features:
        f_out.write(f"--- Feature: {feature} ---\n")
        contingency = pd.crosstab(df['Nuc_type'], df[feature])

        try:
            chi2, p, dof, expected = chi2_contingency(contingency)
            f_out.write(f"Chi-square p-value = {p:.4e}\n")

            contingency = contingency.reindex(ordered_nuctypes)
            overall_prop = df[feature].value_counts(normalize=True)
            row_prop = contingency.div(contingency.sum(axis=1), axis=0)
            fold_enrichment = row_prop.div(overall_prop, axis=1)

            f_out.write("Fold enrichment (row proportion / overall proportion):\n")
            f_out.write(fold_enrichment.round(2).to_string())
            f_out.write("\n\n")

            # Reorder columns if feature in defined orders
            col_order = feature_order_dict.get(feature, fold_enrichment.columns.tolist())
            fold_enrichment = fold_enrichment[col_order]

            # Save heatmap
            plt.figure(figsize=(10, max(3, len(fold_enrichment) * 0.5)))
            sns.heatmap(fold_enrichment, annot=True, cmap='coolwarm', center=1, fmt='.2f',
                        cbar_kws={'label': 'Fold Enrichment'})
            plt.title(f'Fold Enrichment Heatmap: Nuc_type vs {feature}')
            plt.ylabel('Nuc_type')
            plt.xlabel(feature)
            plt.tight_layout()
            plt.savefig(f"enrichment_tables/{feature}_fold_enrichment.png")
            plt.close()  # close the figure to free memory

        except Exception as e:
            f_out.write(f"Could not test due to error: {e}\n\n")


