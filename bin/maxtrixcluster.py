import sys
import os
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------
# 1. Parse command-line argument
# -----------------------------
if len(sys.argv) != 2:
    sys.exit("Usage: python matrixcluster.py <input.cdt>")

input_file = sys.argv[1]
if not os.path.isfile(input_file):
    sys.exit(f"Error: File '{input_file}' not found.")

# -----------------------------
# 2. Load data
# -----------------------------
df = pd.read_csv(input_file, sep="\t", header=None)
df.columns = ['Site', 'Signal1', 'Signal2', 'Signal3', 'Signal4']

# -----------------------------
# 3. Prepare raw and scaled data
# -----------------------------
raw_signals = df[['Signal1', 'Signal2', 'Signal3', 'Signal4']]
scaler = StandardScaler()
signals_scaled = scaler.fit_transform(raw_signals)

# -----------------------------
# 4. KMeans clustering
# -----------------------------
k = 15
kmeans = KMeans(n_clusters=k, random_state=42)
df['Cluster'] = kmeans.fit_predict(signals_scaled)

# -----------------------------
# 5. Save cluster assignment
# -----------------------------
output_tsv = "sites_kmeans_clustered.tsv"
df.to_csv(output_tsv, sep="\t", index=False)
print(f"[✓] Cluster assignments saved to: {output_tsv}")

# -----------------------------
# 6. Compute cluster means (unscaled)
# -----------------------------
cluster_means = df.groupby('Cluster')[['Signal1', 'Signal2', 'Signal3', 'Signal4']].mean()

# -----------------------------
# 7. Plot heatmap of cluster means
# -----------------------------
sns.set(style="white")
g = sns.clustermap(cluster_means, cmap="viridis", method='average', metric='euclidean', figsize=(8, 6))
plt.suptitle(f"Cluster Mean Signal (Unscaled) — KMeans k={k}", y=1.02)

output_png = "kmeans_cluster_means_raw_heatmap.png"
plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.close()
print(f"[✓] Heatmap saved to: {output_png}")
