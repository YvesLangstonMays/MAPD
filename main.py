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

"""
Silhouette score
"""
from sklearn.metrics import silhouette_score

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


"""
Validation via silhouette and gap

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
best_k_sillhouette = K_range[np.argmax(silhouette_scores)]
print(f"optimal k by silhouette: ", {best_k_sillhouette})

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
