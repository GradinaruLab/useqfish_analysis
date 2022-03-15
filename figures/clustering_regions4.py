# %%

# parameters
n_neighbors = 20
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
    "slc30a3",
    "PHP.B8",  # 211229
    "PHP.Astro",
    "CAP-B10",
    "gad2",
    "PHP.V1",
    "cux2",
    "slc17a7",
    "ctgf",
    "sulf2",
    "foxp1",
    "foxp2",
    "calb2",
    "calb1",
    "sst",
    "pvalb",
    "vip",
    "drd1",
    "tac1",
    "lamp5",
    "drd2",
    "gad1",
    "ly6a",
    "pcp4",
    "crym",
    "rorb",
    #  None, 'pcp4', 'crym', 'rorb',
    "rprm",
    "reln",
    "hsd11b1",
    "tpbg",
    "necab1",
    "tnnt1",
    "crh",
    "prkcd",
    "gdf10",
    "lgi2",
    "ppp1r17",
    "gabra6",
    "trh",
    "chrna6",
    "enpp6",
    "th",
    "sncg",
    "chodl",
    "serpinf1",
    "krt73",
    "gpc3",
    "ccnb1",
    None,
    None,
]

gene_list_ordered = [
    "PHP.N",
    "PHP.eB",
    "PHP.B8",  # 211229
    "PHP.Astro",
    "CAP-B10",
    "PHP.V1",
    "slc30a3",
    "gad2",
    "cux2",
    "slc17a7",
    "ctgf",
    "sulf2",
    "foxp1",
    "foxp2",
    "calb2",
    "calb1",
    "sst",
    "pvalb",
    "vip",
    "drd1",
    "tac1",
    "lamp5",
    "drd2",
    "gad1",
    "ly6a",
    "pcp4",
    "crym",
    "rorb",
    #  None, 'pcp4', 'crym', 'rorb',
    "rprm",
    "reln",
    "hsd11b1",
    "tpbg",
    "necab1",
    "tnnt1",
    "crh",
    "prkcd",
    "gdf10",
    "lgi2",
    "ppp1r17",
    "gabra6",
    "trh",
    "chrna6",
    "enpp6",
    "th",
    "sncg",
    "chodl",
    "serpinf1",
    "krt73",
    "gpc3",
    "ccnb1",
]

# %%
from re import L
import scanpy as sc
import scanpy.external as scex
import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np

from glob import glob
import os

from params import *
from matplotlib.pyplot import rc_context, cm

from converting_anndata import *

from sklearn.cluster import AgglomerativeClustering

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
    # sc.pl.heatmap(
    #     adata_sub,
    #     gene_list_ordered[6:],
    #     groupby=daughter_cluster,
    #     cmap=cmap_z,
    #     dendrogram=True,
    #     swap_axes=True,
    #     use_raw=False,
    #     vmax=vmax,
    #     vmin=vmin
    # )

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

            # subind = (adata_sub.obs['cell_id'][0] == cell)
            # break
            # print(f'subind: {subind}')
            # leiden_new.append(n_clusters+1+adata_sub[subind].obs['leiden'][0])
        # print(f'cell: {cell}, {new_cluster_labels}')

    adata_endo.obs[daughter_cluster] = pd.Categorical(new_cluster_labels)

    # clustering adjustment using dendrogram
    # - compare close clusters with two sided Mann-Whitney U test
    # and if not significant, merge two clusters
    sc.tl.dendrogram(
        adata_endo, groupby=daughter_cluster, use_rep="X_pca", linkage_method="ward"
    )


# %%
# ## collect excel files from all positions
path = "./expression_matrices/211229/"
# xlsx2h5ad(path)

# %%
filepath = os.path.join(path, "expression_matrix.h5ad")
adata = sc.read_h5ad(filepath)

print(adata)
print(adata.var_names)
# print(adata.obs['position'])

print(f">>> total cell number: {adata.n_obs}")

# %%

adata_virus = adata[:, adata.var_names.isin(virus_list)]
adata_endo = adata[:, ~adata.var_names.isin(virus_list)]

# data quality check
sc.pp.calculate_qc_metrics(adata_endo, percent_top=None, inplace=True, log1p=False)
sc.pl.violin(
    adata_endo, ["n_genes_by_counts", "total_counts"], jitter=0.4, multi_panel=True
)
sc.pl.scatter(adata_endo, x="total_counts", y="n_genes_by_counts")

print(adata_endo)

# %%
poor_expression_threshold = np.percentile(
    adata_endo.obs["total_counts"], percentile_filter
)
if (poor_expression_threshold == 0) | (np.isnan(poor_expression_threshold)):
    poor_expression_threshold = 1
print(f"poor_expression_threshold: {poor_expression_threshold}")

sns.histplot(adata_endo.obs["total_counts"])

# %%
cell_subset, _ = sc.pp.filter_cells(
    adata_endo, min_counts=poor_expression_threshold, inplace=False
)
cell_subset2, _ = sc.pp.filter_cells(adata_endo, max_counts=4000, inplace=False)
cell_subset = cell_subset & cell_subset2
gene_subset, number_per_gene = sc.pp.filter_genes(
    adata_endo, min_counts=1, inplace=False
)

adata_endo = adata_endo[cell_subset, gene_subset]
adata_virus = adata_virus[cell_subset, :]
print(f">>> total cells passed the filter: {adata_endo.n_obs}")
sns.histplot(adata_endo.obs["total_counts"])

adata_endo.uns["cell_subset"] = cell_subset

# %%
print(adata_endo)
# %%
# print(adata_endo)
adata_endo.raw = adata_endo
adata_endo_norm = sc.pp.normalize_total(adata_endo, copy=True)
adata_endo_log = sc.pp.log1p(adata_endo_norm, copy=True)
adata_endo_scale = sc.pp.scale(adata_endo_log, copy=True)

# %%
fig, axs = plt.subplots(2, 2, figsize=(9, 9))
gene = "calb1"
idx = adata_endo.var_names.get_loc(gene)
# nonzeros = adata_endo.raw.X[:,idx] > 0
nonzeros = adata_endo.raw.X[:, idx] > -1
sns.histplot(adata_endo.raw.X[nonzeros, idx], ax=axs[0, 0])
sns.histplot(adata_endo_norm.X[nonzeros, idx], ax=axs[0, 1])
sns.histplot(adata_endo_log.X[nonzeros, idx], ax=axs[1, 0])
sns.histplot(adata_endo_scale.X[nonzeros, idx], ax=axs[1, 1])

plt.show()


# %%
sc.pp.normalize_total(adata_endo)
# sc.pp.log1p(adata_endo)
sc.pp.scale(adata_endo)


# %%
for idx in range(adata_endo.n_vars):
    sns.scatterplot(
        x=adata_endo.obs["total_counts"], y=adata_endo.raw.X[:, idx], y_jitter=0.1
    )

# %%
for idx in range(adata_endo.n_vars):
    sns.scatterplot(
        x=adata_endo.obs["total_counts"], y=adata_endo.X[:, idx], y_jitter=0.1
    )

# %%
positions = np.asarray(adata_endo.obs["position"])
positions_int = []
for position in positions:
    positions_int.append(int(position[-2:]))
positions = positions_int

position_clusters = []
for position in positions:
    if (position >= 2) & (position <= 6):
        position_clusters.append("cortex")
    elif (position >= 7) & (position <= 11):
        position_clusters.append("striatum")
    elif (position >= 12) & (position <= 17):
        position_clusters.append("thalamus")
    else:
        position_clusters.append("cerebellum")

adata_endo.obs["position_clusters"] = pd.Categorical(position_clusters)
adata_virus.obs["position_clusters"] = pd.Categorical(position_clusters)

with rc_context({"figure.figsize": (15, 15)}):
    sc.pl.heatmap(
        adata_endo,
        # adata_endo.var_names,
        gene_list_ordered[6:],
        groupby="position_clusters",
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=2,
        vmin=-0.5,
    )

# %%
cluster_labels = ["cortex", "striatum", "thalamus", "cerebellum"]
n_clusters = len(cluster_labels)
n_variants = adata_virus.n_vars

transduction_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
transcription_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)

for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(cluster_labels):
        cell_inds = adata_endo.obs["position_clusters"] == cluster
        n_cells_in_cluster = cell_inds.sum()
        transduction_rate[v, c] = (
            np.count_nonzero(adata_virus[cell_inds, virus].X)
            * 100.0
            / n_cells_in_cluster
        )
        nonzeros = (adata_virus[cell_inds, virus].X != 0).ravel()
        tcpr = np.mean(adata_virus[cell_inds, virus].X[nonzeros].ravel())
        if ~np.isnan(tcpr):
            transcription_rate[v, c] = tcpr

# %%
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import cm

blues_zero_white = LinearSegmentedColormap.from_list(
    "", ["white", *cm.Purples(np.arange(255))]
)
reds_zero_white = LinearSegmentedColormap.from_list(
    "", ["white", *cm.Greens(np.arange(255))]
)

fig, axs = plt.subplots(2, 1, figsize=(18, 9))
sns.heatmap(
    data=transduction_rate,
    cmap=blues_zero_white,
    xticklabels=cluster_labels,
    yticklabels=virus_list,
    ax=axs[0],
    linewidths=0.5,
    vmin=0,
    annot=True,
)
sns.heatmap(
    data=transcription_rate,
    cmap=reds_zero_white,
    xticklabels=cluster_labels,
    yticklabels=virus_list,
    ax=axs[1],
    linewidths=0.5,
    vmin=0,
    annot=True,
)
# fig.savefig('./figures/region_broad_virus.svg', format='svg')

# %%
adata_virus.raw = adata_virus
sc.pp.normalize_total(adata_virus)
sc.pp.log1p(adata_virus)
sc.pp.scale(adata_virus)

img = (
    sc.pl.MatrixPlot(
        adata_virus,
        virus_list,
        "position_clusters",
        categories_order=cluster_labels,
        use_raw=False,
        vmax=0.5,
        vmin=-0.5,
        # standard_scale='var'
    )
    .style(cmap=cmap_z)
    .swap_axes()
)
img.savefig("./figures/region_broad_bias.svg", format="svg")

adata_virus = adata_virus.raw.to_adata()


# %%
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
subclustering(
    target_cluster=0, mother_cluster="leiden", daughter_cluster="leiden_sub_0"
)

# %%
subclustering(
    target_cluster=1, mother_cluster="leiden_sub_0", daughter_cluster="leiden_sub_1"
)

# %%
subclustering(
    target_cluster=2, mother_cluster="leiden_sub_1", daughter_cluster="leiden_sub_2"
)


# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
print(adata_virus)

clustering_method = "leiden_sub_2"

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
target_n_clusters = 39
n_merging_branches = 30
# for iter in range(linkage_matrix.shape[0] - target_n_clusters+1):
for iter in range(n_merging_branches):
    cluster1 = int(linkage_matrix[iter, 0])
    cluster2 = int(linkage_matrix[iter, 1])

    print(cluster1, cluster2)
    # inds1 = np.argwhere(cluster_labels==cluster1).ravel()
    # inds2 = np.argwhere(cluster_labels==cluster2).ravel()

    # print(len(inds1), len(inds2))
    # X1 = adata_endo[inds1].raw.X
    # X2 = adata_endo[inds2].raw.X

    # res = mannwhitneyu(X1, X2)

    # print(f'pvalue of cluster {cluster1} and cluster {cluster2}: {res}')
    cluster_labels[np.argwhere(cluster_labels == cluster1).ravel()] = n_clusters
    cluster_labels[np.argwhere(cluster_labels == cluster2).ravel()] = n_clusters
    # print(cluster_labels)

    # if np.any(res.pvalue > p_threshold):
    # if linkage_matrix[iter,2] < 1:
    # print(list(new_labels.values()).index(str(cluster1)))

    cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster1)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster2)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    # new_labels[cinds] = str(n_clusters)
    # new_labels[str(list(new_labels.values()).index(str(cluster2)))] = str(n_clusters)
    # new_labels[str(cluster2)] = str(n_clusters)
    # print(f'clusters merged: {cluster1}+{cluster2} -> {n_clusters}')
    # print(new_labels)

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
# position vs leiden cluster (position-file vs cluster heatmap)
clustering_method = "leiden_new"
positions = list(np.unique(adata_endo.obs["position"]))
cluster_labels = list(np.unique(adata_endo.obs[clustering_method]))
n_pos = len(positions)
n_labels = len(cluster_labels)

position_by_cluster = np.zeros((n_pos, n_labels), dtype=np.uint)

print(positions, cluster_labels)

for cell in range(adata_endo.n_obs):
    pos_ind = positions.index(adata_endo[cell].obs["position"][0])
    cluster_ind = cluster_labels.index(adata_endo[cell].obs[clustering_method][0])
    position_by_cluster[pos_ind, cluster_ind] = (
        position_by_cluster[pos_ind, cluster_ind] + 1
    )

sns.heatmap(
    data=position_by_cluster,
    cmap="viridis",
    xticklabels=cluster_labels,
    yticklabels=positions,
).set(xlabel=clustering_method, ylabel="position")

# %%
total_position_count_per_cluster = np.sum(position_by_cluster, axis=0)

position_by_cluster_norm = (
    position_by_cluster.T / total_position_count_per_cluster[:, None]
)
position_by_cluster_norm = position_by_cluster_norm.T

s = sns.heatmap(
    data=position_by_cluster_norm,
    cmap="viridis",
    xticklabels=cluster_labels,
    yticklabels=positions,
)
s.set(xlabel=clustering_method, ylabel="position")

# %%
s = sns.clustermap(
    data=position_by_cluster_norm,
    cmap="viridis",
    yticklabels=positions,
    row_cluster=False,
    # method='complete',
    # vmax=0.5
)
# s.set(xlabel='leiden', ylabel='position')
# s.fig.savefig('./figures/position_cluster.svg', dpi=600)

# # %%
# leiden_idx_reordered_int = s.dendrogram_col.reordered_ind
# leiden_idx_reordered = [str(idx) for idx in leiden_idx_reordered_int]
# adata_endo.obs['leiden_new2'] = adata_endo.obs['leiden_new'].cat.reorder_categories(leiden_idx_reordered)

# %%
# transduction/transcription rate calculation
cluster_labels = list(np.unique(adata_endo.obs["leiden_new"]))
n_clusters = len(cluster_labels)
n_variants = adata_virus.n_vars

transduction_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
transcription_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)

for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(cluster_labels):
        cell_inds = adata_endo.obs["leiden_new"] == cluster
        n_cell_in_cluster = cell_inds.sum()
        transduction_rate[v, c] = (
            np.count_nonzero(adata_virus[cell_inds, virus].X)
            * 100.0
            / n_cell_in_cluster
        )
        nonzeros = (adata_virus[cell_inds, virus].X != 0).ravel()
        tcpr = np.mean(adata_virus[cell_inds, virus].X[nonzeros].ravel())
        if ~np.isnan(tcpr):
            transcription_rate[v, c] = tcpr


# %%
sc.pl.matrixplot(
    adata_endo,
    gene_list_ordered[6:],
    groupby="leiden_new",
    cmap=cmap_z,
    # dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
)
fig, axs = plt.subplots(2, 1, figsize=(18, 9))
sns.heatmap(
    data=transduction_rate,
    cmap=blues_zero_white,
    xticklabels=list(np.unique(adata_endo.obs["leiden_new"])),
    yticklabels=virus_list,
    ax=axs[0],
    annot=True,
)
sns.heatmap(
    data=transcription_rate,
    cmap=reds_zero_white,
    xticklabels=list(np.unique(adata_endo.obs["leiden_new"])),
    yticklabels=virus_list,
    ax=axs[1],
    annot=True,
)


# %%
#########################
## each region analysis
# striatum

xlsx2h5ad("./expression_matrices/211229/striatum")

# %%

gene_list_striatum = [
    "PHP.N",
    "PHP.eB",
    "PHP.B8",
    "PHP.Astro",
    "CAP-B10",
    "PHP.V1",
    "drd1",
    "drd2",
    "gad2",
    "foxp1",
    "pcp4",
    "tac1",
    "calb1",
    "crym",
    "th",
    "reln",
]
path = "./expression_matrices/211229/striatum/expression_matrix.h5ad"
adata = sc.read_h5ad(path)

print(adata)
print(adata.var_names)
print(f">>> total cell number: {adata.n_obs}")

# %%
adata_virus = adata[:, adata.var_names.isin(virus_list)]
adata_endo = adata[:, adata.var_names.isin(gene_list_striatum[6:])]

sc.pp.calculate_qc_metrics(adata_endo, percent_top=None, inplace=True, log1p=False)
sc.pl.violin(
    adata_endo, ["n_genes_by_counts", "total_counts"], jitter=0.4, multi_panel=True
)
sc.pl.scatter(adata_endo, x="total_counts", y="n_genes_by_counts")

print(adata_endo)

# %%
poor_expression_threshold = np.percentile(
    adata_endo.obs["total_counts"], percentile_filter
)
if (poor_expression_threshold < 1) | (np.isnan(poor_expression_threshold)):
    poor_expression_threshold = 1
print(f"poor_expression_threshold: {poor_expression_threshold}")

sns.histplot(adata_endo.obs["total_counts"])

# %%
cell_subset, _ = sc.pp.filter_cells(
    adata_endo, min_counts=poor_expression_threshold, inplace=False
)
cell_subset2, _ = sc.pp.filter_cells(adata_endo, max_counts=2000, inplace=False)
cell_subset = cell_subset & cell_subset2
gene_subset, number_per_gene = sc.pp.filter_genes(
    adata_endo, min_counts=1, inplace=False
)

adata_endo = adata_endo[cell_subset, gene_subset]
adata_virus = adata_virus[cell_subset, :]
print(f">>> total cells passed the filter: {adata_endo.n_obs}")
sns.histplot(adata_endo.obs["total_counts"])

adata_endo.uns["cell_subset"] = cell_subset

# %%
adata_endo.raw = adata_endo
sc.pp.normalize_total(adata_endo)
# sc.pp.log1p(adata_endo)
sc.pp.scale(adata_endo)

for idx in range(adata_endo.n_vars):
    sns.scatterplot(
        x=adata_endo.obs["total_counts"], y=adata_endo.raw.X[:, idx], y_jitter=0.1
    )

# %%
for idx in range(adata_endo.n_vars):
    sns.scatterplot(
        x=adata_endo.obs["total_counts"], y=adata_endo.X[:, idx], y_jitter=0.1
    )


# %%
sc.pp.pca(adata_endo)
sc.pl.pca_overview(adata_endo)
sc.pp.neighbors(adata_endo, use_rep="X_pca", n_neighbors=n_neighbors)
sc.tl.umap(adata_endo)
sc.tl.tsne(adata_endo, use_rep="X_pca")
sc.tl.leiden(adata_endo, resolution=leiden_resolution)
sc.tl.dendrogram(adata_endo, groupby="leiden", use_rep="X_pca", linkage_method="ward")

# %%
with rc_context({"figure.figsize": (18, 15)}):
    s = sc.pl.heatmap(
        adata_endo,
        gene_list_striatum[6:],
        groupby="leiden",
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=1,
        vmin=-1,
    )
# %%
subclustering(
    target_cluster=0, mother_cluster="leiden", daughter_cluster="leiden_sub_0"
)
with rc_context({"figure.figsize": (18, 15)}):
    s = sc.pl.heatmap(
        adata_endo,
        gene_list_striatum[6:],
        groupby="leiden_sub_0",
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=1,
        vmin=-1,
    )

# %%
subclustering(
    target_cluster=1, mother_cluster="leiden_sub_0", daughter_cluster="leiden_sub_1"
)
with rc_context({"figure.figsize": (18, 15)}):
    s = sc.pl.heatmap(
        adata_endo,
        gene_list_striatum[6:],
        groupby="leiden_sub_1",
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=1,
        vmin=-1,
    )

# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
print(adata_virus)

clustering_method = "leiden_sub_1"

# %%
dendrogram_name = "dendrogram_" + clustering_method
linkage_matrix = adata_endo.uns[dendrogram_name]["linkage"]

cluster_list = [str(idx) for idx in range(adata_endo.obs[clustering_method].nunique())]
adata_endo.obs[clustering_method] = adata_endo.obs[
    clustering_method
].cat.rename_categories(cluster_list)

cluster_labels = np.asarray(adata_endo.obs[clustering_method]).astype(np.uint)
n_clusters = np.unique(cluster_labels).size
new_labels = {str(i): str(i) for i in range(n_clusters)}
p_threshold = 1e-3

sns.lineplot(
    x=np.array(range(linkage_matrix.shape[0])),
    y=linkage_matrix[:, 2],
    marker="o",
)

# %%
arget_n_clusters = 39
n_merging_branches = 14
# for iter in range(linkage_matrix.shape[0] - target_n_clusters+1):
for iter in range(n_merging_branches):
    cluster1 = int(linkage_matrix[iter, 0])
    cluster2 = int(linkage_matrix[iter, 1])

    print(cluster1, cluster2)
    # inds1 = np.argwhere(cluster_labels==cluster1).ravel()
    # inds2 = np.argwhere(cluster_labels==cluster2).ravel()

    # print(len(inds1), len(inds2))
    # X1 = adata_endo[inds1].raw.X
    # X2 = adata_endo[inds2].raw.X

    # res = mannwhitneyu(X1, X2)

    # print(f'pvalue of cluster {cluster1} and cluster {cluster2}: {res}')
    cluster_labels[np.argwhere(cluster_labels == cluster1).ravel()] = n_clusters
    cluster_labels[np.argwhere(cluster_labels == cluster2).ravel()] = n_clusters
    # print(cluster_labels)

    # if np.any(res.pvalue > p_threshold):
    # if linkage_matrix[iter,2] < 1:
    # print(list(new_labels.values()).index(str(cluster1)))

    cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster1)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster2)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    # new_labels[cinds] = str(n_clusters)
    # new_labels[str(list(new_labels.values()).index(str(cluster2)))] = str(n_clusters)
    # new_labels[str(cluster2)] = str(n_clusters)
    # print(f'clusters merged: {cluster1}+{cluster2} -> {n_clusters}')
    # print(new_labels)

    n_clusters = n_clusters + 1

adata_endo.obs["leiden_new"] = (
    adata_endo.obs[clustering_method].map(new_labels).astype("category")
)
# adata_endo.obs['leiden_new'] = adata_endo.obs['leiden']
print(adata_endo.obs["leiden_new"])
sc.tl.dendrogram(
    adata_endo, groupby="leiden_new", use_rep="X_pca", linkage_method="ward"
)

sc.pl.heatmap(
    adata_endo,
    gene_list_striatum[6:],
    groupby="leiden_new",
    cmap=cmap_z,
    dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
)
sc.pl.matrixplot(
    adata_endo,
    gene_list_striatum[6:],
    groupby="leiden_new",
    cmap=cmap_z,
    dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
)
# %%
leiden_idx_reordered_int = [28, 32, 34, 37, 25, 36, 20, 19, 31, 26]
gene_list_striatum = [
    "drd1",
    "tac1",
    "th",
    "foxp1",
    "calb1",
    "drd2",
    "reln",
    "crym",
    "pcp4",
    "gad2",
]
leiden_idx_reordered = [str(idx) for idx in leiden_idx_reordered_int]

# %%
cluster_labels = list(np.unique(adata_endo.obs["leiden_new"]))
# cluster_labels = leiden_idx_reordered
n_clusters = len(cluster_labels)
n_variants = adata_virus.n_vars

transduction_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
transcription_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)

for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(cluster_labels):
        cell_inds = adata_endo.obs["leiden_new"] == cluster
        n_cell_in_cluster = cell_inds.sum()
        transduction_rate[v, c] = (
            np.count_nonzero(adata_virus[cell_inds, virus].X)
            * 100.0
            / n_cell_in_cluster
        )
        nonzeros = (adata_virus[cell_inds, virus].X != 0).ravel()
        tcpr = np.mean(adata_virus[cell_inds, virus].X[nonzeros].ravel())
        if ~np.isnan(tcpr):
            transcription_rate[v, c] = tcpr


# %%
img = (
    sc.pl.MatrixPlot(
        adata_endo,
        gene_list_striatum,
        "leiden_new",
        categories_order=leiden_idx_reordered,
        use_raw=False,
        vmax=vmax,
        vmin=-1,
    )
    .style(cmap="cividis")
    .swap_axes()
)
img.savefig("./figures/region_striatum_endo.svg", format="svg")

fig, axs = plt.subplots(2, 1, figsize=(18, 9))
sns.heatmap(
    data=transduction_rate,
    cmap=blues_zero_white,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[0],
    # annot=True
    linewidths=0.5,
    vmin=0,
)
sns.heatmap(
    data=transcription_rate,
    cmap=reds_zero_white,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[1],
    # annot=True
    linewidths=0.5,
    vmin=0,
)
fig.savefig("./figures/region_striatum_virus.svg", format="svg")

# %%
adata_virus.raw = adata_virus
sc.pp.normalize_total(adata_virus)
sc.pp.log1p(adata_virus)
sc.pp.scale(adata_virus)

bias_striatum = np.zeros((n_variants, n_clusters))
for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(leiden_idx_reordered):
        cell_inds = adata_endo.obs["leiden_new"] == cluster
        bias_striatum[v, c] = np.mean(adata_virus[cell_inds, virus].X.ravel())

img = (
    sc.pl.MatrixPlot(
        adata_virus,
        virus_list,
        groupby="leiden_new",
        categories_order=leiden_idx_reordered,
        use_raw=False,
        vmax=1,
        vmin=-1,
        # standard_scale='var'
    )
    .style(cmap=cmap_z)
    .swap_axes()
)
# img.show()
img.savefig("./figures/region_striatum_bias.svg", format="svg")

adata_virus = adata_virus.raw.to_adata()


# %%
#########################
## each region analysis
# thalamus

xlsx2h5ad("./expression_matrices/211229/thalamus")

# %%
gene_list_thalamus = [
    "PHP.N",
    "PHP.eB",
    "PHP.B8",
    "PHP.Astro",
    "CAP-B10",
    "PHP.V1",
    "rprm",
    "calb1",
    "slc17a7",
    "pvalb",
    "tnnt1",
    "prkcd",
    "calb2",
    "rorb",
    "necab1",
    "pcp4",
]
path = "./expression_matrices/211229/thalamus/expression_matrix.h5ad"
adata = sc.read_h5ad(path)

print(adata)
print(adata.var_names)
print(f">>> total cell number: {adata.n_obs}")

# %%
adata_virus = adata[:, adata.var_names.isin(virus_list)]
adata_endo = adata[:, adata.var_names.isin(gene_list_thalamus[6:])]

sc.pp.calculate_qc_metrics(adata_endo, percent_top=None, inplace=True, log1p=False)
sc.pl.violin(
    adata_endo, ["n_genes_by_counts", "total_counts"], jitter=0.4, multi_panel=True
)
sc.pl.scatter(adata_endo, x="total_counts", y="n_genes_by_counts")

print(adata_endo)

# %%
poor_expression_threshold = np.percentile(
    adata_endo.obs["total_counts"], percentile_filter
)
if (poor_expression_threshold < 1) | (np.isnan(poor_expression_threshold)):
    poor_expression_threshold = 1
print(f"poor_expression_threshold: {poor_expression_threshold}")

sns.histplot(adata_endo.obs["total_counts"])

# %%
cell_subset, _ = sc.pp.filter_cells(
    adata_endo, min_counts=poor_expression_threshold, inplace=False
)
cell_subset2, _ = sc.pp.filter_cells(adata_endo, max_counts=1000, inplace=False)
cell_subset = cell_subset & cell_subset2
gene_subset, number_per_gene = sc.pp.filter_genes(
    adata_endo, min_counts=1, inplace=False
)

adata_endo = adata_endo[cell_subset, gene_subset]
adata_virus = adata_virus[cell_subset, :]
print(f">>> total cells passed the filter: {adata_endo.n_obs}")
sns.histplot(adata_endo.obs["total_counts"])

adata_endo.uns["cell_subset"] = cell_subset

# %%
adata_endo.raw = adata_endo
sc.pp.normalize_total(adata_endo)
# sc.pp.log1p(adata_endo)
sc.pp.scale(adata_endo)

for idx in range(adata_endo.n_vars):
    sns.scatterplot(
        x=adata_endo.obs["total_counts"], y=adata_endo.raw.X[:, idx], y_jitter=0.1
    )

# %%
for idx in range(adata_endo.n_vars):
    sns.scatterplot(
        x=adata_endo.obs["total_counts"], y=adata_endo.X[:, idx], y_jitter=0.1
    )


# %%
sc.pp.pca(adata_endo)
sc.pl.pca_overview(adata_endo)
sc.pp.neighbors(adata_endo, use_rep="X_pca", n_neighbors=n_neighbors)
sc.tl.umap(adata_endo)
sc.tl.tsne(adata_endo, use_rep="X_pca")
sc.tl.leiden(adata_endo, resolution=leiden_resolution)
sc.tl.dendrogram(adata_endo, groupby="leiden", use_rep="X_pca", linkage_method="ward")

# %%
with rc_context({"figure.figsize": (18, 15)}):
    s = sc.pl.heatmap(
        adata_endo,
        gene_list_thalamus[6:],
        groupby="leiden",
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=1,
        vmin=-1,
    )
# %%
subclustering(
    target_cluster=0, mother_cluster="leiden", daughter_cluster="leiden_sub_0"
)
with rc_context({"figure.figsize": (18, 15)}):
    s = sc.pl.heatmap(
        adata_endo,
        gene_list_thalamus[6:],
        groupby="leiden_sub_0",
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=1,
        vmin=-1,
    )

# %%
subclustering(
    target_cluster=1, mother_cluster="leiden_sub_0", daughter_cluster="leiden_sub_1"
)
with rc_context({"figure.figsize": (18, 15)}):
    s = sc.pl.heatmap(
        adata_endo,
        gene_list_thalamus[6:],
        groupby="leiden_sub_1",
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=1,
        vmin=-1,
    )

# %%
subclustering(
    target_cluster=2, mother_cluster="leiden_sub_1", daughter_cluster="leiden_sub_2"
)
with rc_context({"figure.figsize": (18, 15)}):
    s = sc.pl.heatmap(
        adata_endo,
        gene_list_thalamus[6:],
        groupby="leiden_sub_2",
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=1,
        vmin=-1,
    )
# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
print(adata_virus)

clustering_method = "leiden_sub_2"

# %%
dendrogram_name = "dendrogram_" + clustering_method
linkage_matrix = adata_endo.uns[dendrogram_name]["linkage"]

cluster_list = [str(idx) for idx in range(adata_endo.obs[clustering_method].nunique())]
adata_endo.obs[clustering_method] = adata_endo.obs[
    clustering_method
].cat.rename_categories(cluster_list)

cluster_labels = np.asarray(adata_endo.obs[clustering_method]).astype(np.uint)
n_clusters = np.unique(cluster_labels).size
new_labels = {str(i): str(i) for i in range(n_clusters)}
p_threshold = 1e-3

sns.lineplot(
    x=np.array(range(linkage_matrix.shape[0])),
    y=linkage_matrix[:, 2],
    marker="o",
)

# %%
arget_n_clusters = 39
n_merging_branches = 22
# for iter in range(linkage_matrix.shape[0] - target_n_clusters+1):
for iter in range(n_merging_branches):
    cluster1 = int(linkage_matrix[iter, 0])
    cluster2 = int(linkage_matrix[iter, 1])

    print(cluster1, cluster2)
    # inds1 = np.argwhere(cluster_labels==cluster1).ravel()
    # inds2 = np.argwhere(cluster_labels==cluster2).ravel()

    # print(len(inds1), len(inds2))
    # X1 = adata_endo[inds1].raw.X
    # X2 = adata_endo[inds2].raw.X

    # res = mannwhitneyu(X1, X2)

    # print(f'pvalue of cluster {cluster1} and cluster {cluster2}: {res}')
    cluster_labels[np.argwhere(cluster_labels == cluster1).ravel()] = n_clusters
    cluster_labels[np.argwhere(cluster_labels == cluster2).ravel()] = n_clusters
    # print(cluster_labels)

    # if np.any(res.pvalue > p_threshold):
    # if linkage_matrix[iter,2] < 1:
    # print(list(new_labels.values()).index(str(cluster1)))

    cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster1)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster2)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    # new_labels[cinds] = str(n_clusters)
    # new_labels[str(list(new_labels.values()).index(str(cluster2)))] = str(n_clusters)
    # new_labels[str(cluster2)] = str(n_clusters)
    # print(f'clusters merged: {cluster1}+{cluster2} -> {n_clusters}')
    # print(new_labels)

    n_clusters = n_clusters + 1

adata_endo.obs["leiden_new"] = (
    adata_endo.obs[clustering_method].map(new_labels).astype("category")
)
# adata_endo.obs['leiden_new'] = adata_endo.obs['leiden']
print(adata_endo.obs["leiden_new"])
sc.tl.dendrogram(
    adata_endo, groupby="leiden_new", use_rep="X_pca", linkage_method="ward"
)

sc.pl.heatmap(
    adata_endo,
    gene_list_thalamus[6:],
    groupby="leiden_new",
    cmap=cmap_z,
    dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
)
sc.pl.matrixplot(
    adata_endo,
    gene_list_thalamus[6:],
    groupby="leiden_new",
    cmap=cmap_z,
    dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
)


# %%
leiden_idx_reordered_int = [53, 38, 47, 48, 42, 51, 49, 45, 36, 28]
gene_list_thalamus = [
    "prkcd",
    "tnnt1",
    "pvalb",
    "slc17a7",
    "necab1",
    "pcp4",
    "rorb",
    "calb2",
    "calb1",
    "rprm",
]
leiden_idx_reordered = [str(idx) for idx in leiden_idx_reordered_int]
# leiden_idx_reordered = list(np.unique(adata_endo.obs['leiden_new']))

# %%
# cluster_labels = list(np.unique(adata_endo.obs['leiden_new']))
# leiden_idx_reordered = cluster_labels
cluster_labels = leiden_idx_reordered
n_clusters = len(cluster_labels)
n_variants = adata_virus.n_vars

transduction_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
transcription_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)

for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(cluster_labels):
        cell_inds = adata_endo.obs["leiden_new"] == cluster
        n_cell_in_cluster = cell_inds.sum()
        transduction_rate[v, c] = (
            np.count_nonzero(adata_virus[cell_inds, virus].X)
            * 100.0
            / n_cell_in_cluster
        )
        nonzeros = (adata_virus[cell_inds, virus].X != 0).ravel()
        tcpr = np.mean(adata_virus[cell_inds, virus].X[nonzeros].ravel())
        if ~np.isnan(tcpr):
            transcription_rate[v, c] = tcpr


# %%
img = (
    sc.pl.MatrixPlot(
        adata_endo,
        gene_list_thalamus,
        "leiden_new",
        categories_order=leiden_idx_reordered,
        use_raw=False,
        vmax=vmax,
        vmin=-1,
    )
    .style(cmap="cividis")
    .swap_axes()
)
img.savefig("./figures/region_thalamus_endo.svg", format="svg")

fig, axs = plt.subplots(2, 1, figsize=(18, 9))
sns.heatmap(
    data=transduction_rate,
    cmap=blues_zero_white,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[0],
    # annot=True,
    linewidths=0.5,
    vmin=0,
)
sns.heatmap(
    data=transcription_rate,
    cmap=reds_zero_white,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[1],
    # annot=True,
    linewidths=0.5,
    vmin=0,
)
fig.savefig("./figures/region_thalamus_virus.svg", format="svg")

# %%
adata_virus.raw = adata_virus
sc.pp.normalize_total(adata_virus)
sc.pp.log1p(adata_virus)
sc.pp.scale(adata_virus)

bias_thalamus = np.zeros((n_variants, n_clusters))
for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(leiden_idx_reordered):
        cell_inds = adata_endo.obs["leiden_new"] == cluster
        bias_thalamus[v, c] = np.mean(adata_virus[cell_inds, virus].X.ravel())

img = (
    sc.pl.MatrixPlot(
        adata_virus,
        virus_list,
        groupby="leiden_new",
        categories_order=leiden_idx_reordered,
        use_raw=False,
        vmax=1,
        vmin=-1,
        # standard_scale='var'
    )
    .style(cmap=cmap_z)
    .swap_axes()
)
# img.show()
img.savefig("./figures/region_thalamus_bias.svg", format="svg")

adata_virus = adata_virus.raw.to_adata()


# %%
#########################
## each region analysis
# cerebellum

xlsx2h5ad("./expression_matrices/211229/cerebellum")

# %%
gene_list_cerebellum = [
    "PHP.N",
    "PHP.eB",
    "PHP.B8",
    "PHP.Astro",
    "CAP-B10",
    "PHP.V1",
    # 'gad1', 'gad2',
    "prkcd",  #'pvalb',
    "pcp4",
    "calb1",
    "ppp1r17",
    "gdf10",
    "gabra6",
    "calb2",
    "slc17a7",
    "lgi2",
]
path = "./expression_matrices/211229/thalamus/expression_matrix.h5ad"
adata = sc.read_h5ad(path)

print(adata)
print(adata.var_names)
print(f">>> total cell number: {adata.n_obs}")

# %%
adata_virus = adata[:, adata.var_names.isin(virus_list)]
adata_endo = adata[:, adata.var_names.isin(gene_list_cerebellum[6:])]

sc.pp.calculate_qc_metrics(adata_endo, percent_top=None, inplace=True, log1p=False)
sc.pl.violin(
    adata_endo, ["n_genes_by_counts", "total_counts"], jitter=0.4, multi_panel=True
)
sc.pl.scatter(adata_endo, x="total_counts", y="n_genes_by_counts")

print(adata_endo)

# %%
poor_expression_threshold = np.percentile(
    adata_endo.obs["total_counts"], percentile_filter
)
if (poor_expression_threshold < 1) | (np.isnan(poor_expression_threshold)):
    poor_expression_threshold = 1
print(f"poor_expression_threshold: {poor_expression_threshold}")

sns.histplot(adata_endo.obs["total_counts"])

# %%
cell_subset, _ = sc.pp.filter_cells(
    adata_endo, min_counts=poor_expression_threshold, inplace=False
)
cell_subset2, _ = sc.pp.filter_cells(adata_endo, max_counts=1000, inplace=False)
cell_subset = cell_subset & cell_subset2
gene_subset, number_per_gene = sc.pp.filter_genes(
    adata_endo, min_counts=1, inplace=False
)

adata_endo = adata_endo[cell_subset, gene_subset]
adata_virus = adata_virus[cell_subset, :]
print(f">>> total cells passed the filter: {adata_endo.n_obs}")
sns.histplot(adata_endo.obs["total_counts"])

adata_endo.uns["cell_subset"] = cell_subset

# %%
adata_endo.raw = adata_endo
sc.pp.normalize_total(adata_endo)
# sc.pp.log1p(adata_endo)
sc.pp.scale(adata_endo)

for idx in range(adata_endo.n_vars):
    sns.scatterplot(
        x=adata_endo.obs["total_counts"], y=adata_endo.raw.X[:, idx], y_jitter=0.1
    )

# %%
for idx in range(adata_endo.n_vars):
    sns.scatterplot(
        x=adata_endo.obs["total_counts"], y=adata_endo.X[:, idx], y_jitter=0.1
    )


# %%
sc.pp.pca(adata_endo)
sc.pl.pca_overview(adata_endo)
sc.pp.neighbors(adata_endo, use_rep="X_pca", n_neighbors=n_neighbors)
sc.tl.umap(adata_endo)
sc.tl.tsne(adata_endo, use_rep="X_pca")
sc.tl.leiden(adata_endo, resolution=leiden_resolution)
sc.tl.dendrogram(adata_endo, groupby="leiden", use_rep="X_pca", linkage_method="ward")

# %%
with rc_context({"figure.figsize": (18, 15)}):
    s = sc.pl.heatmap(
        adata_endo,
        gene_list_cerebellum[6:],
        groupby="leiden",
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=1,
        vmin=-1,
    )
# %%
subclustering(
    target_cluster=1, mother_cluster="leiden", daughter_cluster="leiden_sub_1"
)
with rc_context({"figure.figsize": (18, 15)}):
    s = sc.pl.heatmap(
        adata_endo,
        gene_list_cerebellum[6:],
        groupby="leiden_sub_1",
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=1,
        vmin=-1,
    )

# # %%
# subclustering(target_cluster=2, mother_cluster='leiden_sub_1', daughter_cluster='leiden_sub_2')
# with rc_context({'figure.figsize': (18, 15)}):
#     s = sc.pl.heatmap(
#         adata_endo,
#         gene_list_cerebellum[6:],
#         groupby='leiden_sub_2',
#         cmap=cmap_z,
#         dendrogram=True,
#         swap_axes=True,
#         use_raw=False,
#         vmax=1,
#         vmin=-1
#     )

# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
print(adata_virus)

clustering_method = "leiden_sub_1"

# %%
dendrogram_name = "dendrogram_" + clustering_method
linkage_matrix = adata_endo.uns[dendrogram_name]["linkage"]

cluster_list = [str(idx) for idx in range(adata_endo.obs[clustering_method].nunique())]
adata_endo.obs[clustering_method] = adata_endo.obs[
    clustering_method
].cat.rename_categories(cluster_list)

cluster_labels = np.asarray(adata_endo.obs[clustering_method]).astype(np.uint)
n_clusters = np.unique(cluster_labels).size
new_labels = {str(i): str(i) for i in range(n_clusters)}
p_threshold = 1e-3

sns.lineplot(
    x=np.array(range(linkage_matrix.shape[0])),
    y=linkage_matrix[:, 2],
    marker="o",
)
# %%
arget_n_clusters = 39
n_merging_branches = 8
# for iter in range(linkage_matrix.shape[0] - target_n_clusters+1):
for iter in range(n_merging_branches):
    cluster1 = int(linkage_matrix[iter, 0])
    cluster2 = int(linkage_matrix[iter, 1])

    # inds1 = np.argwhere(cluster_labels==cluster1).ravel()
    # inds2 = np.argwhere(cluster_labels==cluster2).ravel()

    # print(len(inds1), len(inds2))
    # X1 = adata_endo[inds1].raw.X
    # X2 = adata_endo[inds2].raw.X

    # res = mannwhitneyu(X1, X2)

    # print(f'pvalue of cluster {cluster1} and cluster {cluster2}: {res}')
    cluster_labels[np.argwhere(cluster_labels == cluster1).ravel()] = n_clusters
    cluster_labels[np.argwhere(cluster_labels == cluster2).ravel()] = n_clusters
    # print(cluster_labels)

    # if np.any(res.pvalue > p_threshold):
    # if linkage_matrix[iter,2] < 1:
    # print(list(new_labels.values()).index(str(cluster1)))

    cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster1)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    cinds = [i for i, cluster in new_labels.items() if cluster == str(cluster2)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    # new_labels[cinds] = str(n_clusters)
    # new_labels[str(list(new_labels.values()).index(str(cluster2)))] = str(n_clusters)
    # new_labels[str(cluster2)] = str(n_clusters)
    # print(f'clusters merged: {cluster1}+{cluster2} -> {n_clusters}')
    # print(new_labels)
    print(f"{cluster1} + {cluster2} -> {n_clusters}")
    n_clusters = n_clusters + 1

adata_endo.obs["leiden_new"] = (
    adata_endo.obs[clustering_method].map(new_labels).astype("category")
)
# adata_endo.obs['leiden_new'] = adata_endo.obs['leiden']
print(adata_endo.obs["leiden_new"])
sc.tl.dendrogram(adata_endo, groupby="leiden_new", use_rep="X_pca")

sc.pl.heatmap(
    adata_endo,
    gene_list_cerebellum[6:],
    groupby="leiden_new",
    cmap=cmap_z,
    dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
)
sc.pl.matrixplot(
    adata_endo,
    gene_list_cerebellum[6:],
    groupby="leiden_new",
    cmap=cmap_z,
    dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
)

# %%
leiden_idx_reordered_int = [23, 8, 14, 24, 4, 11, 2, 21, 20]
gene_list_cerebellum = [
    "prkcd",
    "ppp1r17",
    "calb1",
    "pcp4",
    "gdf10",
    "calb2",
    "gabra6",
    "slc17a7",
    "lgi2",
]
# leiden_idx_reordered_int = [
#     38, 7, 6, 37, 40, 39, 9, 11, 10, 32, 1, 35, 31
# ]
# gene_list_cerebellum = [
#     'prkcd', 'pvalb', 'gad1',
#     'pcp4', 'calb1', 'ppp1r17',
#     'gad2', 'gdf10',
#     'calb2', 'gabra6', 'slc17a7', 'lgi2'
# ]
leiden_idx_reordered = [str(idx) for idx in leiden_idx_reordered_int]
# %%
# cluster_labels = list(np.unique(adata_endo.obs['leiden_new']))
cluster_labels = leiden_idx_reordered
n_clusters = len(cluster_labels)
n_variants = adata_virus.n_vars

transduction_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
transcription_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)

for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(cluster_labels):
        cell_inds = adata_endo.obs["leiden_new"] == cluster
        n_cell_in_cluster = cell_inds.sum()
        transduction_rate[v, c] = (
            np.count_nonzero(adata_virus[cell_inds, virus].X)
            * 100.0
            / n_cell_in_cluster
        )
        nonzeros = (adata_virus[cell_inds, virus].X != 0).ravel()
        tcpr = np.mean(adata_virus[cell_inds, virus].X[nonzeros].ravel())
        if ~np.isnan(tcpr):
            transcription_rate[v, c] = tcpr


# %%
img = (
    sc.pl.MatrixPlot(
        adata_endo,
        gene_list_cerebellum,
        "leiden_new",
        categories_order=leiden_idx_reordered,
        use_raw=False,
        vmax=vmax,
        vmin=-1,
    )
    .style(cmap="cividis")
    .swap_axes()
)
img.savefig("./figures/region_cerebellum_endo.svg", format="svg")

fig, axs = plt.subplots(2, 1, figsize=(18, 9))
sns.heatmap(
    data=transduction_rate,
    cmap=blues_zero_white,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[0],
    # annot=True
    linewidths=0.5,
    vmin=0,
)
sns.heatmap(
    data=transcription_rate,
    cmap=reds_zero_white,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[1],
    # annot=True
    linewidths=0.5,
    vmin=0,
)
fig.savefig("./figures/region_cerebellum_virus.svg", format="svg")

# %%
adata_virus.raw = adata_virus
sc.pp.normalize_total(adata_virus)
sc.pp.log1p(adata_virus)
sc.pp.scale(adata_virus)

bias_cerebellum = np.zeros((n_variants, n_clusters))
for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(leiden_idx_reordered):
        cell_inds = adata_endo.obs["leiden_new"] == cluster
        bias_cerebellum[v, c] = np.mean(adata_virus[cell_inds, virus].X.ravel())

img = (
    sc.pl.MatrixPlot(
        adata_virus,
        virus_list,
        groupby="leiden_new",
        categories_order=leiden_idx_reordered,
        use_raw=False,
        vmax=1,
        vmin=-1,
        # standard_scale='var'
    )
    .style(cmap=cmap_z)
    .swap_axes()
)
# img.show()
img.savefig("./figures/region_cerebellum_bias.svg", format="svg")

adata_virus = adata_virus.raw.to_adata()


# %%
