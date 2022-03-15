# %%

# parameters
n_neighbors = 20  # 20
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

gene_list = [  # 220116
    "PHP.N",
    "PHP.eB",
    "slc30a3",
    "PHP.B8",
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
    "gad1",
    "tac1",
    "lamp5",
    "crym",
    "rorb",
    "rprm",
    "reln",
    "hsd11b1",
    "tpbg",
    "necab1",
    "pcp4",
    "crh",
    "th",
    "trh",
    "lgi2",
    "ppp1r17",
    "krt73",
    "sncg",
    "chrna6",
    None,
    None,
]

gene_list_ordered = [
    "PHP.eB",
    "CAP-B10",
    "PHP.N",
    "PHP.Astro",
    "PHP.V1",
    "PHP.B8",
    "slc30a3",
    "slc17a7",
    "cux2",
    "calb1",
    "lamp5",
    "foxp1",
    "rorb",
    "sulf2",
    "th",
    "pcp4",
    "rprm",
    "crym",
    "foxp2",
    "necab1",
    "hsd11b1",
    "tpbg",
    "ctgf",
    "trh",
    "chrna6",
    "ppp1r17",
    "gad1",
    "gad2",
    "pvalb",
    "reln",
    "sst",
    "tac1",
    "calb2",
    "crh",
    "vip",
    "sncg",
    "krt73",
    "lgi2",
]

import os
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
    sc.pp.neighbors(adata_sub)
    sc.tl.leiden(adata_sub, key_added=daughter_cluster)
    sc.tl.dendrogram(adata_sub, groupby=daughter_cluster)
    sc.pl.heatmap(
        adata_sub,
        gene_list_ordered[6:],
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

            # subind = (adata_sub.obs['cell_id'][0] == cell)
            # break
            # print(f'subind: {subind}')
            # leiden_new.append(n_clusters+1+adata_sub[subind].obs['leiden'][0])
        # print(f'cell: {cell}, {new_cluster_labels}')

    adata_endo.obs[daughter_cluster] = pd.Categorical(new_cluster_labels)

    # clustering adjustment using dendrogram
    # - compare close clusters with two sided Mann-Whitney U test
    # and if not significant, merge two clusters
    sc.tl.dendrogram(adata_endo, groupby=daughter_cluster, use_rep="X_pca")
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


from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import cm

blues_zero_white = LinearSegmentedColormap.from_list(
    "", ["white", *cm.Blues(np.arange(255))]
)
reds_zero_white = LinearSegmentedColormap.from_list(
    "", ["white", *cm.Reds(np.arange(255))]
)

## collect excel files from all positions
# path = './expression_matrices/210828/220117_analyzed'
path = "./expression_matrices/220116/expression_matrix_stitched.h5ad"
adata = sc.read_h5ad(path)

print(adata)
print(adata.var_names)
# print(adata.obs['position'])

print(f">>> total cell number: {adata.n_obs}")

# %%

adata_virus = adata[~adata.obs["overlap_in_stitched"], adata.var_names.isin(virus_list)]
adata_endo = adata[~adata.obs["overlap_in_stitched"], ~adata.var_names.isin(virus_list)]

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
print(adata_endo)
# %%
# print(adata_endo)
adata_endo.raw = adata_endo
# %%
# adata_endo.X = adata_endo.X / adata_endo.X.sum(axis=0)[None,:] * adata_endo.n_vars
# adata_endo.X = adata_endo.X / adata_endo.X.sum(axis=1)[:,None] * adata_endo.n_obs
# adata_endo.X = adata_endo.X / adata_endo.X.sum(axis=0)[None,:]
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
sc.pp.pca(adata_endo)
sc.pl.pca_overview(adata_endo)
sc.pp.neighbors(adata_endo, use_rep="X_pca", n_neighbors=n_neighbors)
sc.tl.umap(adata_endo)
sc.tl.tsne(adata_endo, use_rep="X_pca")
sc.tl.leiden(adata_endo, resolution=leiden_resolution)
sc.tl.dendrogram(
    adata_endo, groupby="leiden", use_rep="X_pca", linkage_method="complete"
)

# %%
with rc_context({"figure.figsize": (18, 15)}):
    s = sc.pl.heatmap(
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

# # %%
subclustering(
    target_cluster=2, mother_cluster="leiden_sub_1", daughter_cluster="leiden_sub_2"
)

# # %%
subclustering(
    target_cluster=3, mother_cluster="leiden_sub_2", daughter_cluster="leiden_sub_3"
)
# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
print(adata_virus)

clustering_method = "leiden_sub_3"
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
    # sc.pl.heatmap(
    #     adata_endo,
    #     gene_list_ordered[6:],
    #     groupby=clustering_method,
    #     cmap=cmap_raw,
    #     dendrogram=True,
    #     swap_axes=True,
    #     use_raw=True,
    #     standard_scale='var'
    # )
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
# # position vs leiden cluster (position-file vs cluster heatmap)
# positions = np.asarray(adata_endo.obs['position'])
# positions_int = []
# for position in positions:
#     positions_int.append(int(position[-2:]))
# positions = positions_int
# labels = np.asarray(adata_endo.obs[clustering_method]).astype(np.int)
# n_pos = np.amax(positions)+1
# n_labels = np.amax(labels)+1

# position_by_cluster = np.zeros((n_pos, n_labels), dtype=np.int)

# for pos, label in zip(positions, labels):
#     position_by_cluster[pos, label] = position_by_cluster[pos,label] + 1

# # active_position = np.sum(position_by_cluster, axis=1) > 0
# active_position = [
#     0, 11, 12, 23,
#     1, 10, 13, 22,
#     2, 9, 14, 21,
#     3, 8, 15, 20,
#     4, 7, 16, 19,
#     5, 6, 17, 18,
# ]
# position_by_cluster = position_by_cluster[np.array(active_position),:]
# # active_position = np.argwhere(active_position).flatten()
# # print(active_position)
# s = sns.heatmap(data=position_by_cluster, cmap='viridis', yticklabels=active_position)
# s.set(xlabel=clustering_method, ylabel='position')

# %%
# total_position_count_per_cluster = np.sum(position_by_cluster, axis=0)

# position_by_cluster_norm = position_by_cluster.T/total_position_count_per_cluster[:,None]
# position_by_cluster_norm = position_by_cluster_norm.T

# s = sns.heatmap(
#     data=position_by_cluster_norm,
#     cmap='viridis',
#     yticklabels=active_position
# )
# s.set(xlabel=clustering_method, ylabel='position')

# s = sns.clustermap(
#     data=position_by_cluster_norm,
#     cmap='viridis',
#     yticklabels=active_position,
#     row_cluster=False,
#     method='ward',
#     # vmax=0.5
# )
# # s.set(xlabel='leiden', ylabel='position')
# # s.fig.savefig('./figures/position_cluster.svg', dpi=600)


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
n_merging_branches = 44
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
# cluster ordering

# leiden_idx_reordered_int = [
#     9, 10, 6, 0, 18, 17, 15, 13, 5, 12, 16, 2, 3, 14, 7, 1, 4, 8, 11  # with neighbor=50
# ]
# leiden_idx_reordered_int = [
#     10, 8, 0, 14, 17, 2, 15, 3, 12, 11, 23, 28, 27, 1, 16, 20, 6, 13, 5, 4, 9, 7, 18    #with neighbor=20, dendrogram merging
# ]
leiden_idx_reordered_int = [
    96,  # foxp1
    97,  # lgi2
    95,  # foxp2
    115,  # chrna6/hsd11b1
    6,  # trh
    103,  # necab1                                         # ppp1r17
    108,  # slc17a7
    24,
    110,
    109,
    79,
    114,
    32,
    86,  # L2/3/4
    105,  # L4
    62,
    104,  # L5
    112,
    101,
    85,
    82,  # L6
    51,  # Sp
    3,
    107,
    40,
    29,
    71,
    91,
]

gene_list_ordered = [
    "PHP.eB",
    "CAP-B10",
    "PHP.N",
    "PHP.Astro",
    "PHP.V1",
    "PHP.B8",
    "foxp1",
    "lgi2",
    "foxp2",
    "chrna6",
    "trh",
    "necab1",
    "slc17a7",
    "lamp5",
    "slc30a3",
    "cux2",
    "calb1",
    "ppp1r17",
    "rorb",
    "calb2",
    "crh",
    "hsd11b1",
    "pcp4",
    "tpbg",
    "crym",
    "rprm",
    "ctgf",
    "sulf2",
    "th",
    "krt73",
    "gad1",
    "gad2",
    "pvalb",
    "tac1",
    "reln",
    "sst",
    "vip",
    "sncg",
]


leiden_idx_reordered = [str(idx) for idx in leiden_idx_reordered_int]
# leiden_idx_reordered = list(np.unique(adata_endo.obs['leiden_new']))
# adata_endo.obs['leiden_new'] = adata_endo.obs['leiden_new'].cat.reorder_categories(leiden_idx_reordered)

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
# cluster_labels = list(np.unique(adata_endo.obs['leiden_new']))
n_clusters = len(cluster_labels)
n_variants = adata_virus.n_vars

transduction_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
transcription_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
mean_count = np.zeros((2, n_clusters))

for v, virus in enumerate(virus_list):
    for c, cluster in enumerate(leiden_idx_reordered):
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
        if v == 0:
            mean_count[0, c] = np.mean(adata_endo[cell_inds].raw.X.ravel())
            mean_count[1, c] = np.sum(adata_endo[cell_inds].raw.X.ravel())

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
sns.heatmap(data=mean_count, xticklabels=leiden_idx_reordered, annot=True, ax=axs[2])

# %%
sc.pp.normalize_total(adata_virus)
sc.pp.log1p(adata_virus)
sc.pp.scale(adata_virus)

# %%
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
img.savefig("./figures/layer_endo.svg", format="svg")

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
img.savefig("./figures/layer_bias.svg", format="svg")

fig, axs = plt.subplots(2, 1, figsize=(18, 9))
sns.heatmap(
    data=transduction_rate,
    # cmap=blues_zero_white,
    cmap="viridis",
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[0],
    linewidths=0.5,
    vmin=0,
)
sns.heatmap(
    data=transcription_rate,
    cmap=reds_zero_white,
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list,
    ax=axs[1],
    linewidths=0.5,
    vmin=0,
)
fig.savefig("./figures/layer_virus.svg", format="svg")


# %%
# visualization
clustering_method = "leiden_new"
n_clusters = np.unique(list(adata_endo.obs[clustering_method])).size

cell_labels_stitched = zarr.load(
    "./expression_matrices/220116/cell_labels_stitched.zarr"
)
cell_cmap = get_cmap("rainbow", n_clusters)
cell_colors = [rgb2hex(cell_cmap(i)) for i in range(n_clusters)]
cell_coords = []
for cluster in leiden_idx_reordered:
    print(f"cluster: {cluster}")
    cell_locations_str = np.asarray(
        adata_endo[(adata_endo.obs[clustering_method] == cluster)].obs["cell_location"]
    )
    print(len(cell_locations_str))
    cell_locations = []
    for cell_loc in cell_locations_str:
        cell_loc_float = [float(coord) for coord in re.findall(r"[\d]*[.]+", cell_loc)]
        cell_locations.append(cell_loc_float)

    print(np.array(cell_locations).shape)
    # print(cell_locations)
    cell_coords.append(np.array(cell_locations))

viewer = napari.Viewer()

viewer.add_labels(cell_labels_stitched, name="cell_labels", multiscale=False)

for c, cell_coord in enumerate(cell_coords):
    viewer.add_points(
        cell_coord,
        face_color=cell_colors[c],
        edge_color=cell_colors[c],
        size=100,
        n_dimensional=True,
        name=f"cluster {leiden_idx_reordered[c]}",
        blending="additive",
        visible=False,
    )


# # %%
# fig, axs = plt.subplots(6, 4, figsize=(18,24))
# for c, cell_coord in enumerate(cell_coords):
#     axs[c].scatter(cell_coords[:,1], cell_coord[:,0])
# %%
ind = 26
viewer.add_points(
    cell_coords[ind],
    face_color="#AA4499",
    edge_color="#AA4499",
    size=100,
    n_dimensional=True,
    name=f"cluster {leiden_idx_reordered[ind]}",
    blending="additive",
)
# %%

# %%
