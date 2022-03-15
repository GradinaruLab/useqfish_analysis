# %%
# # parameters
n_neighbors = 30  # 20
vmax = 3
vmin = -3
cmap_z = "coolwarm"  #'RdBu_r'
cmap_raw = "viridis"
ifdoublet = False
percentile_filter = 30
min_counts = 3
poor_cluster_percentile_threshold = 20
poor_cluster_zstd_threshold = 10
leiden_resolution = 0.5

gene_list = [
    "PHP.N",
    "PHP.eB",
    "PHP.B8",  # 210828
    "PHP.Astro",
    "CAP-B10",
    "PHP.V1",
    "gad1",
    "slc17a7",
    "gja1",  #'mbp',
    "cldn5",
    "hexb",
    "msr1",
    "acta2",
    "sst",
    "pvalb",
    "vip",
]

gene_list_ordered = [
    "PHP.eB",
    "CAP-B10",
    "PHP.N",
    "PHP.Astro",
    "PHP.V1",
    "PHP.B8",
    "slc17a7",
    "gad1",
    "pvalb",
    "sst",
    "vip",
    "gja1",
    "hexb",
    "msr1",
    "cldn5",
    "acta2",
]

import re
from glob import glob

import matplotlib.pyplot as plt
import napari
import numpy as np

# %%
import scanpy as sc
import scanpy.external as scex
import seaborn as sns
import zarr
from matplotlib.cm import get_cmap
from matplotlib.colors import rgb2hex
from matplotlib.pyplot import rc_context, savefig
from sklearn.cluster import AgglomerativeClustering

from converting_anndata import *
from params import *

from warnings import filterwarnings

filterwarnings("ignore")

sc.settings.verbosity = 3


def anndata_copy_attributes(src_anndata, dtn_anndata):
    dtn_anndata.obs = src_anndata.obs
    dtn_anndata.uns = src_anndata.uns
    dtn_anndata.obsm = src_anndata.obsm
    dtn_anndata.obsp = src_anndata.obsp


def subclustering(target_cluster, mother_cluster, daughter_cluster):
    # subclustering
    adata_sub = adata_endo[adata_endo.obs[mother_cluster] == str(target_cluster), :]
    print(adata_sub)
    sc.pp.pca(adata_sub)
    sc.pp.neighbors(adata_sub, use_rep="X_pca")
    sc.tl.leiden(adata_sub, key_added=daughter_cluster, resolution=leiden_resolution)
    sc.tl.dendrogram(
        adata_sub, use_rep="X_pca", groupby=daughter_cluster, linkage_method="ward"
    )
    sc.pl.heatmap(
        adata_sub,
        gene_list_ordered[7:],
        groupby=daughter_cluster,
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin,
    )

    cluster_labels = np.asarray(adata_endo.obs[mother_cluster]).astype(np.uint)

    n_clusters = np.unique(cluster_labels).size
    new_cluster_labels = []

    for cell in range(adata_endo.n_obs):
        # print(adata_endo[cell].obs['leiden'])
        if adata_endo[cell].obs[mother_cluster][0] != str(target_cluster):
            new_cluster_labels.append(adata_endo[cell].obs[mother_cluster][0])
        else:
            cell_id = adata_endo[cell].obs[mother_cluster].index
            subcluster_label = adata_sub.obs[daughter_cluster].loc[cell_id][0]
            new_cluster_labels.append(str(n_clusters + 1 + int(subcluster_label)))

    adata_endo.obs[daughter_cluster] = pd.Categorical(new_cluster_labels)

    # clustering adjustment using dendrogram
    # - compare close clusters with two sided Mann-Whitney U test
    # and if not significant, merge two clusters
    sc.tl.dendrogram(
        adata_endo, groupby=daughter_cluster, use_rep="X_pca", linkage_method="ward"
    )
    sc.pl.heatmap(
        adata_endo,
        gene_list_ordered[6:],
        groupby=daughter_cluster,
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin,
    )
    sc.pl.matrixplot(
        adata_endo,
        gene_list_ordered[6:],
        groupby=daughter_cluster,
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin,
    )


# %%
path = "./expression_matrices/major"
# path = './expression_matrices/210430'
filepath = os.path.join(path, "*.h5ad")
filenames = sorted(glob(filepath), key=os.path.basename)
adatas = [sc.read_h5ad(filename) for filename in filenames]
adata = sc.concat(adatas, join="inner")
adata.obs_names_make_unique()

print(adata)
print(adata.var_names)
# print(adata.obs['position'])

print(f">>> total cell number: {adata.n_obs}")

# %%
# endogenous and viral gene separation
adata_virus = adata[:, adata.var_names.isin(virus_list)]
adata_endo = adata[:, adata.var_names.isin(gene_list[6:])]

# data quality check
sc.pp.calculate_qc_metrics(adata_endo, percent_top=None, inplace=True, log1p=False)
sc.pl.violin(
    adata_endo, ["n_genes_by_counts", "total_counts"], jitter=0.4, multi_panel=True
)
sc.pl.scatter(adata_endo, x="total_counts", y="n_genes_by_counts")

# %%
# data filtering
poor_expression_threshold = np.percentile(adata_endo.obs["total_counts"], 20)
if poor_expression_threshold == 0:
    poor_expression_threshold = 1
print(f"poor_expression_threshold: {poor_expression_threshold}")

sns.histplot(adata_endo.obs["total_counts"])

# %%
cell_subset, _ = sc.pp.filter_cells(
    adata_endo, min_counts=poor_expression_threshold, inplace=False
)
cell_subset2, _ = sc.pp.filter_cells(adata_endo, max_counts=1e3, inplace=False)
cell_subset = cell_subset & cell_subset2
gene_subset, number_per_gene = sc.pp.filter_genes(
    adata_endo, min_counts=3, inplace=False
)

adata_endo = adata_endo[cell_subset, gene_subset]
adata_virus = adata_virus[cell_subset, :]
print(f">>> total cells passed the filter: {adata_endo.n_obs}")

# %%
# normalization
# adata_endo.X = adata_endo.X / adata_endo.X.sum(axis=0)[None,:] * adata_endo.n_vars
# adata_endo.X = adata_endo.X / adata_endo.X.sum(axis=1)[:,None] * adata_endo.n_obs

adata_endo.raw = adata_endo
adata_endo_norm = sc.pp.normalize_total(adata_endo, copy=True)
adata_endo_log = sc.pp.log1p(adata_endo_norm, copy=True)
adata_endo_scale = sc.pp.scale(adata_endo_log, copy=True)

# %%
fig, axs = plt.subplots(3, 2, figsize=(9, 9))
idx = 1
# nonzeros = adata_endo.raw.X[:,idx] > 0
nonzeros = adata_endo.raw.X[:, idx] > -1
sns.scatterplot(
    x=adata_endo.obs["total_counts"],
    y=adata_endo.raw.X[:, idx],
    y_jitter=0.1,
    ax=axs[0, 0],
)
sns.histplot(adata_endo.raw.X[nonzeros, idx], ax=axs[0, 1])
sns.histplot(adata_endo_norm.X[nonzeros, idx], ax=axs[1, 0])
sns.histplot(adata_endo_log.X[nonzeros, idx], ax=axs[1, 1])
sns.histplot(adata_endo_scale.X[nonzeros, idx], ax=axs[2, 0])
sns.scatterplot(
    x=adata_endo.obs["total_counts"],
    y=adata_endo_scale.X[:, idx],
    y_jitter=0.1,
    ax=axs[2, 1],
)

plt.show()

# %%
sc.pp.normalize_total(adata_endo)
# sc.pp.log1p(adata_endo)
sc.pp.scale(adata_endo)

# %%
for idx in range(adata_endo.n_vars):
    # idx = 11
    sns.scatterplot(
        x=adata_endo.obs["total_counts"], y=adata_endo.raw.X[:, idx], y_jitter=0.1
    )

# %%
for idx in range(adata_endo.n_vars):
    sns.scatterplot(
        x=adata_endo.obs["total_counts"], y=adata_endo.X[:, idx], y_jitter=0.1
    )

# %%
# data clustering
sc.pp.pca(adata_endo)
sc.pl.pca_overview(adata_endo)
sc.pp.neighbors(adata_endo, use_rep="X_pca", n_neighbors=n_neighbors)
sc.tl.umap(adata_endo)
sc.tl.tsne(adata_endo, use_rep="X_pca")
sc.tl.leiden(adata_endo, resolution=leiden_resolution)
sc.tl.dendrogram(adata_endo, groupby="leiden", use_rep="X_pca", linkage_method="ward")

# %%
sc.pl.heatmap(
    adata_endo,
    gene_list_ordered[6:],
    groupby="leiden",
    cmap=cmap_z,
    dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
)
sc.pl.matrixplot(
    adata_endo,
    gene_list_ordered[6:],
    groupby="leiden",
    cmap=cmap_z,
    dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
)
# %%
# subclustering(target_cluster=7, mother_cluster='leiden', daughter_cluster='leiden_sub_7')

# %%
# subclustering(target_cluster=19, mother_cluster='leiden_sub_7', daughter_cluster='leiden_sub_19')
# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
print(adata_virus)

clustering_method = "leiden"

# %%
sc.pl.umap(
    adata_virus,
    color=["PHP.eB", "CAP-B10", "PHP.N", "PHP.Astro", "PHP.B8", "PHP.V1", "leiden"],
    frameon=True,
    ncols=3,
    vmax="p99",
)
sc.pl.tsne(
    adata_virus,
    color=["PHP.eB", "CAP-B10", "PHP.N", "PHP.Astro", "PHP.B8", "PHP.V1", "leiden"],
    frameon=True,
    ncols=3,
    vmax="p99",
)
sc.pl.tsne(adata_endo, color=["slc17a7", "gad1", "leiden"], frameon=True, vmax="p99")

# %%
with rc_context({"figure.figsize": (20, 20)}):
    sc.pl.heatmap(
        adata_endo,
        # adata_endo.var_names,
        gene_list_ordered[6:],
        groupby=clustering_method,
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin,
    )
    sc.pl.matrixplot(
        adata_endo,
        gene_list_ordered[6:],
        groupby=clustering_method,
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin,
    )

# %%
dendrogram_name = "dendrogram_" + clustering_method
linkage_matrix = adata_endo.uns[dendrogram_name]["linkage"]
print(linkage_matrix)

cluster_list = [str(idx) for idx in range(adata_endo.obs[clustering_method].nunique())]
adata_endo.obs[clustering_method] = adata_endo.obs[
    clustering_method
].cat.rename_categories(cluster_list)

cluster_labels = np.asarray(adata_endo.obs[clustering_method]).astype(np.uint)
n_clusters = np.unique(cluster_labels).size
new_labels = {str(i): str(i) for i in range(n_clusters)}
p_threshold = 1e-3

# %%
sns.lineplot(
    x=np.array(range(linkage_matrix.shape[0])),
    y=linkage_matrix[:, 2],
    marker="o",
)
# %%
n_merging_branches = 4
# for iter in range(linkage_matrix.shape[0] - target_n_clusters+1):
for iter in range(n_merging_branches):
    cluster1 = int(linkage_matrix[iter, 0])
    cluster2 = int(linkage_matrix[iter, 1])

    print(cluster1, cluster2)

    cluster_labels[np.argwhere(cluster_labels == cluster1).ravel()] = n_clusters
    cluster_labels[np.argwhere(cluster_labels == cluster2).ravel()] = n_clusters

    cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster1)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster2)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    n_clusters = n_clusters + 1


cluster1 = 6
cluster2 = 17
cluster_labels[np.argwhere(cluster_labels == cluster1).ravel()] = n_clusters
cluster_labels[np.argwhere(cluster_labels == cluster2).ravel()] = n_clusters
cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster1)]
for cind in cinds:
    new_labels[cind] = str(n_clusters)
cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster2)]
for cind in cinds:
    new_labels[cind] = str(n_clusters)
n_clusters = n_clusters + 1

adata_endo.obs["leiden_new"] = (
    adata_endo.obs[clustering_method].map(new_labels).astype("category")
)
# adata_endo.obs['leiden_new'] = adata_endo.obs['leiden']
print(adata_endo.obs["leiden_new"])
sc.tl.dendrogram(adata_endo, groupby="leiden_new", use_rep="X_pca")

sc.pl.matrixplot(
    adata_endo,
    gene_list_ordered[6:],
    groupby="leiden_new",
    cmap=cmap_z,
    dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
)

# %%
leiden_idx_reordered_int = [
    # 17, 18, 14, 15, 8, 16, 6, 13, 12, 11
    # 23, 6, 22, 18, 8, 17, 19, 14, 7, 20
    # 23, 24, 18, 8, 17, 19, 14, 7, 20
    # 16, 6, 17, 14, 8, 15, 7, 11, 12, 10
    16,
    18,
    14,
    8,
    15,
    7,
    11,
    12,
    10,
]

# poor_idx = ['9']
# cells_outlier_cluster = adata_endo.obs['leiden_new'].isin(poor_idx)
# adata_endo = adata_endo[~cells_outlier_cluster,:]
# adata_virus = adata_virus[~cells_outlier_cluster,:]

leiden_idx_reordered = [str(idx) for idx in leiden_idx_reordered_int]
adata_endo.obs["leiden_new"] = adata_endo.obs["leiden_new"].cat.reorder_categories(
    leiden_idx_reordered
)

img = (
    sc.pl.MatrixPlot(
        adata_endo,
        gene_list_ordered[6:],
        "leiden_new",
        categories_order=leiden_idx_reordered,
        use_raw=False,
        vmax=vmax,
        vmin=vmin,
    )
    .style(cmap=cmap_z)
    .swap_axes()
)
img.show()
# %%
# transduction/transcription rate calculation
cluster_labels = list(np.unique(adata_endo.obs["leiden_new"]))
n_clusters = len(cluster_labels)
n_variants = adata_virus.n_vars

transduction_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
transcription_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
# enrichment = np.zeros((n_variants, n_clusters))
# transduction_bias = np.zeros((n_variants, n_clusters))


for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(leiden_idx_reordered):
        cell_inds = adata_endo.obs["leiden_new"] == cluster
        n_cell_in_cluster = cell_inds.sum()
        transduction_rate[v, c] = (
            np.count_nonzero(adata_virus[cell_inds, virus].X)
            * 100.0
            / n_cell_in_cluster
        )
        # enrichment[v,c] = np.mean(np.log1p(adata_virus[cell_inds, virus].X.ravel()))
        nonzeros = (adata_virus[cell_inds, virus].X != 0).ravel()
        tcpr = np.mean(adata_virus[cell_inds, virus].X[nonzeros].ravel())
        if ~np.isnan(tcpr):
            transcription_rate[v, c] = tcpr
        # transduction_bias[v,c] = np.count_nonzero(adata_virus[cell_inds, virus].X)

# transduction_bias = transduction_bias/np.sum(transduction_bias, axis=1)[:,None]
# delta_bias =

from matplotlib.pyplot import cm
from matplotlib.colors import LinearSegmentedColormap

blues_zero_white = LinearSegmentedColormap.from_list(
    "", ["white", *cm.Purples(np.arange(255))]
)
reds_zero_white = LinearSegmentedColormap.from_list(
    "", ["white", *cm.Greens(np.arange(255))]
)

fig, axs = plt.subplots(3, 1, figsize=(18, 9))
sns.heatmap(
    data=transduction_rate,
    cmap=blues_zero_white,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[0],
    annot=True,
)
sns.heatmap(
    data=transcription_rate,
    cmap=reds_zero_white,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[1],
    annot=True,
)

# %%
# adata_virus_norm = sc.pp.normalize_total(adata_virus, copy=True)
# sc.pp.log1p(adata_virus_norm)
adata_virus_norm = adata_virus.copy()
adata_virus_norm.X = adata_virus_norm.X / adata_endo.obs["total_counts"][:, None]

sc.pp.normalize_total(adata_virus)
sc.pp.log1p(adata_virus)
sc.pp.scale(adata_virus)

# sc.pp.normalize_total(adata_virus)
# sc.pp.log1p(adata_virus)
# sc.pp.scale(adata_virus)

# %%
enrichment = np.zeros((n_variants, n_clusters))
bias = np.zeros((n_variants, n_clusters))
for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(leiden_idx_reordered):
        cell_inds = adata_endo.obs["leiden_new"] == cluster
        bias[v, c] = np.mean(adata_virus[cell_inds, virus].X.ravel())
        enrichment[v, c] = np.mean(adata_virus_norm[cell_inds, virus].X.ravel())

# %%
fig, axs = plt.subplots(2, 1, figsize=(18, 9))
sns.heatmap(
    data=bias,
    cmap=cmap_z,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    vmax=0.5,
    vmin=-0.5,
    ax=axs[0],
)
sns.heatmap(
    data=enrichment,
    cmap="viridis",
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    annot=True,
    ax=axs[1],
)

# %%
img = (
    sc.pl.MatrixPlot(
        adata_endo,
        gene_list_ordered[6:],
        "leiden_new",
        categories_order=leiden_idx_reordered,
        use_raw=False,
        vmax=vmax,
        vmin=-1,
    )
    .style(cmap="cividis")
    .swap_axes()
)
img.savefig("./figures/major_endo.svg", format="svg")

img = (
    sc.pl.MatrixPlot(
        adata_virus,
        virus_list,
        "leiden_new",
        categories_order=leiden_idx_reordered,
        use_raw=False,
        vmax=0.5,
        vmin=-0.5,
    )
    .style(cmap=cmap_z)
    .swap_axes()
)
img.savefig("./figures/major_bias.svg", format="svg")

fig, axs = plt.subplots(3, 1, figsize=(18, 18))
sns.heatmap(
    data=transduction_rate,
    cmap=blues_zero_white,
    # cmap='viridis',
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[0],
    linewidths=0.5,
    vmin=0,
    # annot=True
)
sns.heatmap(
    data=transcription_rate,
    cmap=reds_zero_white,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[1],
    linewidths=0.5,
    vmin=0,
    # annot=True
)
sns.heatmap(
    data=enrichment,
    cmap="viridis",
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[2],
    linewidths=0.5,
    vmin=0,
)
fig.savefig("./figures/major_virus.svg", format="svg")

# %%
