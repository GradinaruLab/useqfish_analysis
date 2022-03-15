# %%

# parameters
n_neighbors = 50
vmax = 3
vmin = -1
cmap_z = "coolwarm"
cmap_raw = "viridis"
ifdoublet = False
poor_cluster_percentile_threshold = 20
poor_cluster_zstd_threshold = 10

# %%
import scanpy as sc
import scanpy.external as scex
import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np

from glob import glob
import os

from params import *
from matplotlib.pyplot import rc_context, savefig

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


## collect excel files from all positions
# path = './expression_matrices/210828/220117_analyzed'
path = "./expression_matrices/211229"
filepath = os.path.join(path, "*.h5ad")
filenames = sorted(glob(filepath), key=os.path.basename)
print(filenames)

if len(filenames) == 0:
    xlsx2h5ad(path)

# %%
adatas = [sc.read_h5ad(filename) for filename in filenames]

if len(adatas) == 1:
    adata = adatas[0]
else:
    # for adata in adatas: print(adata.var_names)
    adata = sc.concat(adatas, join="inner")
    adata.obs_names_make_unique()

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

# %%
print(adata_endo)
# %%
poor_expression_threshold = np.percentile(adata_endo.obs["total_counts"], 50)
if (poor_expression_threshold == 0) | (np.isnan(poor_expression_threshold)):
    poor_expression_threshold = 1
print(f"poor_expression_threshold: {poor_expression_threshold}")

sns.histplot(adata_endo.obs["total_counts"])

# %%
cell_subset, _ = sc.pp.filter_cells(
    adata_endo, min_counts=poor_expression_threshold, inplace=False
)
cell_subset2, _ = sc.pp.filter_cells(adata_endo, max_counts=1500, inplace=False)
cell_subset = cell_subset & cell_subset2
gene_subset, number_per_gene = sc.pp.filter_genes(
    adata_endo, min_counts=1, inplace=False
)

adata_endo = adata_endo[cell_subset, gene_subset]
adata_virus = adata_virus[cell_subset, :]
print(f">>> total cells passed the filter: {adata_endo.n_obs}")

# %%
print(adata_endo)
# %%
# print(adata_endo)
adata_endo.raw = adata_endo
adata_endo_norm = sc.pp.normalize_total(adata_endo, target_sum=1e4, copy=True)
adata_endo_log = sc.pp.log1p(adata_endo_norm, copy=True)
adata_endo_scale = sc.pp.scale(adata_endo_log, copy=True)

# %%
fig, axs = plt.subplots(2, 2, figsize=(9, 9))
idx = 5
# nonzeros = adata_endo.raw.X[:,idx] > 0
nonzeros = adata_endo.raw.X[:, idx] > -1
sns.histplot(adata_endo.raw.X[nonzeros, idx], ax=axs[0, 0])
sns.histplot(adata_endo_norm.X[nonzeros, idx], ax=axs[0, 1])
sns.histplot(adata_endo_log.X[nonzeros, idx], ax=axs[1, 0])
sns.histplot(adata_endo_scale.X[nonzeros, idx], ax=axs[1, 1])

plt.show()

# %%
sc.pp.normalize_total(adata_endo, target_sum=1e4)
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
sc.pp.neighbors(adata_endo, use_rep="X_pca", n_neighbors=n_neighbors)
sc.tl.umap(adata_endo)
sc.tl.tsne(adata_endo, use_rep="X_pca")
sc.tl.leiden(adata_endo, resolution=leiden_resolution)
sc.tl.dendrogram(adata_endo, groupby="leiden", use_rep="X_pca")


# %%
# poor clustering adjustment based on std(z-score) across the genes
# cluster_labels = np.asarray(adata_endo.obs['leiden']).astype(np.uint)
# n_clusters = cluster_labels.max()+1
# n_clusters = n_clusters.astype(np.uint)
# stdevs = np.zeros((n_clusters,))
# for cluster in range(n_clusters):
#     inds = np.argwhere(cluster_labels==cluster).ravel()
#     mean_z = np.mean(adata_endo[inds].X, axis=0)
#     stdevs[cluster] = np.std(mean_z)
# threshold = np.percentile(stdevs, poor_cluster_zstd_threshold)
# # threshold = .2
# #
# new_labels = {str(i):str(i) for i in range(n_clusters.astype(np.uint))}
# for cluster, stdev in enumerate(stdevs):
#     if stdev < threshold:
#         new_labels[str(cluster)] = 'poor'

# print(stdevs)
# print(threshold)
# print(new_labels)
# adata_endo.obs['leiden_new'] = adata_endo.obs['leiden'].map(new_labels)


# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
print(adata_virus)


# %%

sc.pl.umap(
    adata_virus,
    color=["PHP.eB", "CAP-B10", "PHP.N", "PHP.Astro", "PHP.B8", "PHP.V1", "leiden"],
    s=50,
    frameon=True,
    ncols=3,
    vmax="p99",
)


sc.pl.tsne(
    adata_virus,
    color=["PHP.eB", "CAP-B10", "PHP.N", "PHP.Astro", "PHP.B8", "PHP.V1", "leiden"],
    s=50,
    frameon=True,
    ncols=3,
    vmax="p99",
)


clustering_method = "leiden"

# %%
with rc_context({"figure.figsize": (15, 15)}):
    sc.pl.heatmap(
        adata_endo,
        adata_endo.var_names,
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
        adata_endo.var_names,
        groupby=clustering_method,
        cmap=cmap_z,
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin,
    )

# %%
cluster_labels = np.asarray(adata_endo.obs[clustering_method]).astype(np.uint)
print(cluster_labels)
n_clusters = np.unique(cluster_labels).size
print(n_clusters)
n_variants = adata_virus.n_vars

transduction_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
transcription_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)

for i, virus in enumerate(virus_list):
    for j, cluster in enumerate(np.unique(cluster_labels)):
        cluster = int(cluster)
        cell_inds = np.argwhere(cluster_labels == cluster).ravel()
        transduction_rate[i, j] = (
            np.count_nonzero(adata_virus[cell_inds, virus].X) * 100.0 / cell_inds.size
        )
        nonzeros = (adata_virus[cell_inds, virus].X != 0).ravel()
        tcpr = np.mean(adata_virus[cell_inds, virus].X[nonzeros].ravel())
        if ~np.isnan(tcpr):
            transcription_rate[i, j] = tcpr


fig, axs = plt.subplots(1, 2, figsize=(18, 9))
sns.heatmap(
    data=transduction_rate,
    cmap="Blues",
    xticklabels=np.unique(cluster_labels),
    yticklabels=virus_list,
    ax=axs[0],
)
sns.heatmap(
    data=transcription_rate,
    cmap="Reds",
    xticklabels=np.unique(cluster_labels),
    yticklabels=virus_list,
    ax=axs[1],
)

# #######################################
# # %%
# with rc_context({'figure.figsize': (15, 15)}):
#     sc.pl.matrixplot(
#         adata_endo,
#         # adata_endo.var_names,
#         gene_list_ordered[6:],
#         'leiden',
#         swap_axes=True,
#         dendrogram=True,
#         use_raw=False,
#         cmap=cmap_z,
#         vmax=vmax,
#         vmin=vmin
#     )

#     sc.pl.matrixplot(
#         adata_virus,
#         adata_virus.var_names,
#         'leiden',
#         cmap='Blues',
#         standard_scale='obs',
#         colorbar_title='column log scaled\nexpression',
#         swap_axes=True,
#         dendrogram=True,
#     )
#     sc.pl.matrixplot(
#         adata_virus,
#         adata_virus.var_names,
#         'leiden',
#         cmap='Reds',
#         standard_scale='var',
#         colorbar_title='row log scaled\nexpression',
#         swap_axes=True,
#         dendrogram=True,
#     )

# # graph after poor cluster filtering
# inds = adata_endo.obs['leiden_new'] != 'poor'
# adata_endo = adata_endo[inds, :]
# adata_virus = adata_virus[inds, :]
# sc.tl.dendrogram(adata_endo, groupby='leiden_new')
# anndata_copy_attributes(adata_endo, adata_virus)

# with rc_context({'figure.figsize': (15, 15)}):
#     sc.pl.heatmap(
#         adata_endo,
#         adata_endo.var_names,
#         groupby='leiden_new',
#         cmap=cmap_z,
#         dendrogram=True,
#         swap_axes=True,
#         use_raw=False,
#         vmax=vmax,
#         vmin=vmin
#         # log=True
#     )

# with rc_context({'figure.figsize': (15, 15)}):
#     sc.pl.matrixplot(
#         adata_endo,
#         # adata_endo.var_names,
#         gene_list_ordered[6:],
#         'leiden_new',
#         swap_axes=True,
#         dendrogram=True,
#         use_raw=False,
#         vmax=vmax,
#         vmin=vmin,
#         cmap=cmap_z
#     )
#     sc.pl.matrixplot(
#         adata_endo,
#         # adata_endo.var_names,
#         gene_list_ordered[6:],
#         'leiden_new',
#         swap_axes=True,
#         dendrogram=True,
#         cmap=cmap_raw
#     )
#     sc.pl.stacked_violin(
#         adata_endo,
#         gene_list_ordered[6:],
#         'leiden_new',
#         dendrogram=True,
#         swap_axes=True
#     )

#     sc.pl.dotplot(
#         adata_virus,
#         adata_virus.var_names,
#         'leiden_new',
#         cmap='Blues',
#         standard_scale='obs',
#         colorbar_title='column scaled\nexpression',
#         swap_axes=True,
#         dendrogram=True,

#     )
#     sc.pl.dotplot(
#         adata_virus,
#         adata_virus.var_names,
#         'leiden_new',
#         cmap='Reds',
#         standard_scale='var',
#         colorbar_title='row scaled\nexpression',
#         swap_axes=True,
#         dendrogram=True,
#     )
#     sc.pl.matrixplot(
#         adata_virus,
#         adata_virus.var_names,
#         'leiden_new',
#         cmap='Reds',
#         standard_scale='var',
#         colorbar_title='row scaled\nexpression',
#         swap_axes=True,
#         dendrogram=True,
#     )


# %%
# position vs leiden cluster (position-file vs cluster heatmap)
positions = np.asarray(adata_endo.obs["position"])
positions_int = []
for position in positions:
    positions_int.append(int(position[-2:]))
positions = positions_int
labels = np.asarray(adata_endo.obs["leiden"]).astype(np.int)
n_pos = np.amax(positions) + 1
n_labels = np.amax(labels) + 1

position_by_cluster = np.zeros((n_pos, n_labels), dtype=np.int)

for pos, label in zip(positions, labels):
    position_by_cluster[pos, label] = position_by_cluster[pos, label] + 1

active_position = np.sum(position_by_cluster, axis=1) > 0
position_by_cluster = position_by_cluster[active_position, :]
active_position = np.argwhere(active_position).flatten()
print(active_position)
s = sns.heatmap(data=position_by_cluster, cmap="viridis", yticklabels=active_position)
s.set(xlabel="leiden", ylabel="position")

# %%
total_position_count_per_cluster = np.sum(position_by_cluster, axis=0)

position_by_cluster_norm = (
    position_by_cluster.T / total_position_count_per_cluster[:, None]
)
position_by_cluster_norm = position_by_cluster_norm.T

s = sns.heatmap(
    data=position_by_cluster_norm, cmap="viridis", yticklabels=active_position
)
s.set(xlabel="leiden", ylabel="position")

# %%
s = sns.clustermap(
    data=position_by_cluster_norm,
    cmap="viridis",
    yticklabels=active_position,
    row_cluster=False,
    method="ward",
    # vmax=0.5
)
s.fig.savefig("./figures/position_cluster.svg", dpi=600)

# %%
leiden_idx_reordered = s.dendrogram_col.reordered_ind
leiden_idx_reordered = [str(idx) for idx in leiden_idx_reordered]
with rc_context({"figure.figsize": (15, 15)}):
    img = sc.pl.MatrixPlot(
        adata_endo,
        # adata_endo.var_names,
        gene_list_ordered[6:],
        "leiden",
        categories_order=leiden_idx_reordered,
        use_raw=False,
        cmap="bwr",
    ).swap_axes()
    img.show()
    img.savefig("./figures/region_endo.svg", format="svg")

    img = sc.pl.MatrixPlot(
        adata_endo,
        gene_list_ordered[6:],
        "leiden",
        categories_order=leiden_idx_reordered,
        cmap="viridis",
    ).swap_axes()
    img.show()
    img.savefig("./figures/region_endo_raw.svg", format="svg")

    img = (
        sc.pl.MatrixPlot(
            adata_virus,
            ["PHP.eB", "CAP-B10", "PHP.N", "PHP.Astro", "PHP.V1", "PHP.B8"],
            "leiden",
            categories_order=leiden_idx_reordered,
            standard_scale="obs",
        )
        .style(cmap="Blues")
        .swap_axes()
    )
    img.show()
    img.savefig("./figures/region_virus_obs.svg", format="svg")

    img = (
        sc.pl.MatrixPlot(
            adata_virus,
            ["PHP.eB", "CAP-B10", "PHP.N", "PHP.Astro", "PHP.V1", "PHP.B8"],
            "leiden",
            categories_order=leiden_idx_reordered,
            standard_scale="var",
        )
        .style(cmap="Reds")
        .swap_axes()
    )
    img.show()
    img.savefig("./figures/region_virus_var.svg", format="svg")
# %%
# sc.tl.rank_genes_groups(adata_endo, groupby='leiden', method='wilcoxon')
# sc.pl.rank_genes_groups_matrixplot(adata_endo, n_genes=1, use_raw=False, vmin=-2, vmax=2, cmap='bwr')
# sc.pl.rank_genes_groups_heatmap(adata_endo, n_genes=1, use_raw=False, vmin=-2, vmax=2, cmap='bwr')


# %%

## cortex - major cell types
# img = sc.pl.MatrixPlot(
#     adata_endo,
#     ['slc17a7', 'gad1', 'pvalb', 'sst', 'vip', 'gja1', 'mbp', 'msr1', 'hexb', 'cldn5', 'acta2'],
#     'leiden',
#     categories_order=['0', '3', '4', '1', '9', '2', '7', '8', '6','5']
# ).swap_axes()
# img.savefig('./figures/cortex_endo.svg', format='svg')

# img = sc.pl.MatrixPlot(
#     adata_virus,
#     ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8'],
#     'leiden',
#     categories_order=['0', '3', '4', '1', '9', '2', '7', '8', '6','5'],
#     standard_scale='obs',
# ).style(cmap='Blues').swap_axes()
# img.savefig('./figures/cortex_virus_obs.svg', format='svg')

# img = sc.pl.MatrixPlot(
#     adata_virus,
#     ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8'],
#     'leiden',
#     categories_order=['0', '3', '4', '1', '9', '2', '7', '8','6','5'],
#     standard_scale='var',
# ).style(cmap='Reds').swap_axes()
# img.savefig('./figures/cortex_virus_var.svg', format='svg')

## cortex - subtypes
# img = sc.pl.MatrixPlot(
#     adata_endo,
#     adata_endo.var_names,
#     'leiden',
#     categories_order=['5', '12', '9', '18', '21', '0', '4', '1', '19', '10', '20', '22', '14', '17', '15', '8','3', '6', '24', '7', '13', '16', '11', '2', '23'],
#     # swap_axes=True,
#     # dendrogram=True
# ).swap_axes()
# img.savefig('./figures/endo.svg', format='svg')

# img = sc.pl.MatrixPlot(
#     adata_virus,
#     ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8'],
#     'leiden',
#     categories_order=['5', '12', '9', '18', '21', '0', '4', '1', '19', '10', '20', '22', '14', '17', '15', '8','3', '6', '24', '7', '13', '16', '11', '2', '23'],
#     standard_scale='obs',
# ).style(cmap='Blues').swap_axes()
# img.savefig('./figures/virus_obs.svg', format='svg')

# img = sc.pl.MatrixPlot(
#     adata_virus,
#     ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8'],
#     'leiden',
#     categories_order=['5', '12', '9', '18', '21', '0','4', '1', '19', '10', '20', '22', '14', '17', '15', '8','3', '6', '24', '7', '13', '16', '11', '2', '23'],
#     standard_scale='var',
# ).style(cmap='Reds').swap_axes()
# img.savefig('./figures/virus_var.svg', format='svg')

## striatum
# img = sc.pl.MatrixPlot(
#     adata_endo,
#     ['tac1', 'tpbg', 'reln', 'crym', 'gad1', 'gad2', 'sst', 'pcp4'],
#     'leiden',
#     categories_order=['2', '3', '4', '7', '8', '0', '1', '9', '5', '6'],
# ).swap_axes()
# img.savefig('./figures/striatum_endo.svg', format='svg')

# img = sc.pl.MatrixPlot(
#     adata_virus,
#     ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8'],
#     'leiden',
#     categories_order=['2', '3', '4', '7', '8', '0', '1', '9', '5', '6'],
#     standard_scale='obs',
# ).style(cmap='Blues').swap_axes()
# img.savefig('./figures/striatum_virus_obs.svg', format='svg')

# img = sc.pl.MatrixPlot(
#     adata_virus,
#     ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8'],
#     'leiden',
#     categories_order=['2', '3', '4', '7', '8', '0', '1', '9', '5', '6'],
#     standard_scale='var',
# ).style(cmap='Reds').swap_axes()
# img.savefig('./figures/striatum_virus_var.svg', format='svg')

## midbrain
# img = sc.pl.MatrixPlot(
#     adata_endo,
#     ['slc6a4', 'gad1', 'th', 'sst'],
#     'leiden',
#     categories_order=['0', '1', '2', '3', '4'],
# ).swap_axes()
# img.savefig('./figures/mid_endo.svg', format='svg')

# img = sc.pl.MatrixPlot(
#     adata_virus,
#     ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8'],
#     'leiden',
#     categories_order=['0', '1', '2', '3', '4'],
#     standard_scale='obs',
# ).style(cmap='Blues').swap_axes()
# img.savefig('./figures/mid_virus_obs.svg', format='svg')

# img = sc.pl.MatrixPlot(
#     adata_virus,
#     ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8'],
#     'leiden',
#     categories_order=['0', '1', '2', '3', '4'],
#     standard_scale='var',
# ).style(cmap='Reds').swap_axes()
# img.savefig('./figures/mid_virus_var.svg', format='svg')
# %%
