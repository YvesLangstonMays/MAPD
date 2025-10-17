# pandas for dataframe manipulations
import pandas as pd

# numpy for numeric operations and specialized arrays
import numpy as np

# plotting library for figures
import matplotlib.pyplot as plt


# pathlib for file management and filesystem navigation
from pathlib import Path

"""
sklearn is a standard ML toolkit
principal component analysis reduces the dimensions of the data linearly
down to only the principal components that explain most of the variation
"""
from sklearn.decomposition import PCA

from tqdm import tqdm

from scipy.stats import entropy

"""
Silhouette score
"""
from sklearn.metrics import silhouette_score

"""
Since our gap gave an optimal k of 7, and the silhouette gave an optimal
k of 2, there may be a hierarchical structure to our data.

These two libraries will help elucidate the hierachy embedded in our data.
"""

from sklearn.mixture import GaussianMixture
import hdbscan

"""
t-distributed stochastic neighbor embedding is nonlinear dimensionality reduction for
visualization
preserves local neighborhoods, so rather than preserving absolute distances,
it preserves the probabilities of two points being neighbors in the original
dimensional space
"""
from sklearn.manifold import TSNE

""" KMeans is a clustering algorithm that assigns each point to one of k centroids by minimizing
squared euclidean distance within the clusters """
from sklearn.cluster import KMeans

results_dir = Path("Results")
results_dir.mkdir(exist_ok=True)

silhouette_summary = {}

# set global k val
k_val = 3

df = pd.read_csv("OlsenData_TableS6.csv")


time_columns = ['"0" EGF', '"1" EGF', '"5" EGF', '"10" EGF', '"20" EGF']

"""
For clustering and distance metrics, missing values would invalidate our calculations,
so we convert our values to numeric, coerce any possible errors to NaN, drop those NaN values
from the dataframe and finally, ensure that our data is of type float32.
"""

X = (
    df[time_columns]
    .apply(pd.to_numeric, errors="coerce")
    .dropna()
    .to_numpy(dtype=np.float32)
)

"""
calculate z score for each peptide, row wise. 
Rows are designated by the axis=1, dimensions are retained with argument keepdims=True
Each peptide's mean is subtracted, and the result is divided by its standard deviation

Using Z score normalized the magnitude differences between peptides so clustering emphasizes
shape of peptide time course instead of the absolute intensity.
Using Z score sets each peptide to a mean of 0 and std of 1 across its 5 time points.
"""

Xz = (X - X.mean(axis=1, keepdims=True)) / (X.std(axis=1, keepdims=True))

"""
Here, principal component analyssi is performed to find orthogonal directions that capture
the largest variance. Projecting into the first two components gives a 2D summary of the main
dynamic patterns

n_components argument ensures we are computing the first two principal components
random_state is used by habit, to make any stochastic pieces deterministic
"""
pca = PCA(n_components=2, random_state=42)
Z = pca.fit_transform(Xz)


plt.figure(figsize=(6, 5))
plt.scatter(Z[:, 0], Z[:, 1], s=8, alpha=0.6)
plt.title(
    f"PCA (PC1={pca.explained_variance_ratio_[0]*100:.1f}%, PC2={pca.explained_variance_ratio_[1]*100:.1f}%)"
)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.savefig(results_dir / f"pca_variance_k_{k_val}.png", dpi=300)
plt.close()

"""
tsne is used here for visualization, and doesn't replace PCA for analysis.

Perplexity represents how many neighbors each point cares about
random_state sets the seed
init 
learning rate is the step size for gradient updates. auto lets sklearn pick a rate based on
data size
"""
tsne = TSNE(
    n_components=2, perplexity=30, random_state=42, init="pca", learning_rate="auto"
)

U = tsne.fit_transform(Z)

"""
Simple k means clustering.. Partitions data int ok clusters by minimizing the 
within cluster squared euclidean distance to centroids

here we start with 3 clusters, run k means with 10 different centroid seeds and pick the best
random_state ensures reproducible labels
"""

kmeans = KMeans(n_clusters=k_val, n_init=10, random_state=42)
labels = kmeans.fit_predict(Z)

plt.figure(figsize=(6, 5))
plt.scatter(U[:, 0], U[:, 1], c=labels, cmap="viridis", s=10)
plt.title("t-SNE of EGF response clusters")
plt.xlabel("t-SNE 1")
plt.ylabel("t-SNE 2")
plt.colorbar(label="Cluster")
plt.savefig(results_dir / f"tsne_clusters_k_{k_val}.png", dpi=300)
plt.close()
"""
X_df includes each column as a time point

"""
minutes = [0, 1, 5, 10, 20]
X_df = pd.DataFrame(Xz, columns=minutes)
X_df["cluster"] = labels
plt.figure(figsize=(7, 5))
for c in sorted(X_df["cluster"].unique()):
    mu = X_df[X_df.cluster == c].iloc[:, :-1].mean()
    sd = X_df[X_df.cluster == c].iloc[:, :-1].std()
    plt.plot(minutes, mu, marker="o", label=f"Cluster {c}")
    plt.fill_between(minutes, mu - sd, mu + sd, alpha=0.2)
plt.xlabel("Time (min)")
plt.ylabel("Z-scored phospho-signal")
plt.title("Average EGF response per cluster")
plt.legend()
plt.savefig(results_dir / f"mean_timecourses_k_{k_val}.png", dpi=300)
plt.close()

################### Silhouette Scores ##################
"""
Validation via silhouette and gap for kmeans

"""

# initialize list
silhouette_scores = []

# define k range for cluster testing
K_range = range(2, 10)

# compute score and add to list
for k in K_range:
    kmeans = KMeans(n_clusters=k, n_init=10, random_state=42)
    labels = kmeans.fit_predict(Z)
    sil = silhouette_score(Z, labels)
    silhouette_scores.append(sil)

plt.figure(figsize=(6, 4))
plt.plot(K_range, silhouette_scores, marker="o")
plt.title("Silhouette Score vs. # of Clusters")
plt.xlabel("# of clusters (k)")
plt.ylabel("Mean Silhouette Score")
plt.grid(True)
plt.tight_layout()
plt.savefig(results_dir / f"silhouette_scores_k_{k_val}.png", dpi=300)
plt.close()

# print best k
best_kmeans_sillhouette = K_range[np.argmax(silhouette_scores)]
silhouette_summary["KMeans_best_k"] = best_kmeans_sillhouette
silhouette_summary["KMeans_silhouette"] = max(silhouette_scores)

"""

Increasing number of clusters reduces the within cluster sum of squared distances, referred tp as
intertia. If inertia is used as a metric, the pitfall of increasing clusters
is inevitable. 

The gap statistic aims to fix this by measuring how much better the clustering is than
what would be expected from random noise.

If the data has real clusters, the inertia will be much lower than random data with no
structure. If adding more cluseters doesnt imrpove much, beyond what random noise would,
the k needs to be lowered.

A gap is the difference betwwen the expected log(Wk) from a reference dataset
and the observed log(Wk) from the actual dataset. 

Both Wk values represent within cluster disperson, applied to the actual and random
datasets, respectively

The optimal number of clusters, k is the smallest k such that the Gap(k) >= Gap(k + 1) - s(k+1)
where s(k+1) is the standard deviation of the simulated gap values that adjusts for simulation 
variability in the reference datasets

Here, we pick the k with the largest gap value since that represents a large difference
between our actual data and random noise.
"""


def gap_stat(data, refs=10, max_k=10):
    gaps = np.zeros(max_k - 1)
    for k in tqdm(range(1, max_k), desc="computing gap stat"):
        km = KMeans(n_clusters=k_val, n_init=10, random_state=42)
        km.fit(data)
        disp = np.log(km.inertia_)

        ref_disps = np.zeros(refs)
        for i in range(refs):
            random_ref = np.random.random_sample(size=data.shape)
            km_ref = KMeans(n_clusters=k_val, n_init=10, random_state=42)
            km_ref.fit(random_ref)
            ref_disps[i] = np.log(km_ref.inertia_)
        # computes the gaps for each k
        gaps[k - 1] = np.mean(ref_disps) - disp
    return gaps


gaps = gap_stat(Z, refs=10, max_k=10)
plt.figure(figsize=(6, 4))
plt.plot(range(1, 10), gaps, marker="o")
plt.title("Gap Stat vs. # of Clusters")
plt.xlabel("# of clusters (k)")
plt.ylabel("Gap Value")
plt.grid(True)
plt.tight_layout()
plt.savefig(results_dir / f"gap_stat_k_{k_val}.png", dpi=300)
plt.close()

best_k_gap = np.argmax(gaps) + 1
print(f"optimal k by gap: {best_k_gap}")

silhouette_summary["Gap_best_k"] = best_k_gap
silhouette_summary["Gap_values"] = gaps.tolist()


############## DB SCAN & GMM ###############
"""
GMM reveals soft clusters (overlapping subgroups).
Each point has a membership probability to all clusters

"""

gmm_silhouette_scores = []
for k in K_range:
    gmm = GaussianMixture(n_components=k, covariance_type="full", random_state=42)
    gmm_labels = gmm.fit_predict(Z)
    sil = silhouette_score(Z, gmm_labels)
    gmm_silhouette_scores.append(sil)

plt.figure(figsize=(6, 4))
plt.plot(K_range, gmm_silhouette_scores, marker="o", color="crimson")
plt.title("GMM  Score vs. # of Components")
plt.xlabel("# of components (k)")
plt.ylabel("GMM Score")
plt.grid(True)
plt.tight_layout()
plt.savefig(results_dir / f"gmm_scores_k_{k_val}.png", dpi=300)
plt.close()

best_k_gmm = K_range[np.argmax(gmm_silhouette_scores)]
silhouette_summary["GMM_best_k"] = best_k_gmm
silhouette_summary["GMM_silhouette"] = max(gmm_silhouette_scores)

# visualize GMM clustering on tSNE space
gmm_best = GaussianMixture(
    n_components=best_k_gmm, covariance_type="full", random_state=42
)
gmm_labels = gmm_best.fit_predict(Z)

plt.figure(figsize=(6, 5))
plt.scatter(U[:, 0], U[:, 1], c=gmm_labels, cmap="plasma", s=10)
plt.title("tSNE of GMM Clusters")
plt.xlabel("tSNE 1")
plt.ylabel("tSNE 2")
plt.colorbar(label="Cluster")
plt.savefig(results_dir / f"tsne_gmm_clusters_k_{k_val}.png", dpi=300)
plt.close()


"""
HDBSCAN revals hard clusters and outliers.

"""

# min_cluster_size defines the smallest size of dense regions
hdb = hdbscan.HDBSCAN(min_cluster_size=20, min_samples=5, cluster_selection_epsilon=0.0)
hdb_labels = hdb.fit_predict(Z)

# compute silhouette only if multiple clusters exist
if len(set(hdb_labels)) > 1:
    hdb_sil = silhouette_score(Z, hdb_labels)
    silhouette_summary["HDBSCAN_silhouette"] = hdb_sil
else:
    silhouette_summary["HDBSCAN_silhouette"] = None

plt.figure(figsize=(6, 5))
plt.scatter(U[:, 0], U[:, 1], c=hdb_labels, cmap="Spectral", s=10)
plt.title("t-SNE of HDBSCAN Clusters")
plt.xlabel("t-SNE 1")
plt.ylabel("t-SNE 2")
plt.colorbar(label="Cluster")
plt.savefig(results_dir / f"tsne_hdbscan_clusters_k_{k_val}.png", dpi=300)
plt.close()

print("\n\n==== Validation Scores Summary  ====")
for method, score in silhouette_summary.items():
    print(f"{method:20s}: {score}")


print("\n\n")
"""
Olsen et al. use of fuzzy was correct, since there are no distinct boundaries between the data.

Exploring memberships as experimental variables
"""

memberships = df.filter(like="membership").apply(pd.to_numeric, errors="coerce")
memberships.sum(axis=1).describe()
memberships.max(axis=1).hist(bins=50)
plt.title("Dist of max membership per peptide")
plt.xlabel("max membership")
plt.ylabel("Freq")
plt.show()

df["dominant_cluster"] = memberships.idxmax(axis=1)
df["dominant_strength"] = memberships.max(axis=1)


for c in memberships.columns:
    mask = df["dominant_cluster"] == c
    mu = df.loc[mask, time_columns].mean()
    plt.plot([0, 1, 5, 10, 20], mu, label=f"{c} (n={mask.sum()})")
plt.legend()
plt.title("avg egf response by dom cluster")
plt.xlabel("Time (min)")
plt.ylabel("Phospho ratio")
plt.show()

"""
'fuzziness index' = entropy per peptide
entropy closer to 0 means the peptide is tightly assigned to one cluster
entropy closer to 1 means the peptide shares membership broadly and there may possible signaling cross-talk
"""
df["entropy"] = memberships.apply(lambda x: entropy(x, base=6), axis=1)

entropy_by_protein = df.groupby("Accession")["entropy"].mean()
print("Entropy by protein \n")
print(entropy_by_protein)

"""
Next: mapping to known pathways
taking the dominant cluster or top 2 memberships per peptide
map them to proteins using uniprot or phosphositeplus
do enrichment analysis, reactome or STRING, per cluster

this should show which clusters correspond to Early EGFR autophosphorylation events, mapk/erk cascades,
cytosekeltal remodeling, nuclear transcriptional feedback

then I need to visualize it somehow.. maybe a membership heatmap, peptides vs clusters, colored by strength
entropy histogram, network overlay where i can map cluster memberships onto a string protein interaction
graph
"""

top_membership = entropy_by_protein.max()
top_membership_index = entropy_by_protein.idxmax()

print(f"\nTop member: ", {top_membership}, "Top member ID: ", {top_membership_index})
# find the sequence corresponding to the top member

ipi_ids = ["IPI00306280"]

peptides = df.loc[df["Accession"] == ipi_id, "Phosphopeptide sequence"].unique()
print(f"Peptides for {ipi_id}:")
for pep in peptides:
    print(pep)
