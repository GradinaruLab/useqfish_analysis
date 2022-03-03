# %%

# parameters
n_neighbors = 50
vmax = 3
vmin = -3
cmap_z = 'coolwarm' #'RdBu_r'
cmap_raw = 'viridis'
ifdoublet = False
percentile_filter = 30
min_counts = 3
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
from matplotlib.pyplot import rc_context, cm

from converting_anndata import *

from sklearn.cluster import AgglomerativeClustering

from warnings import filterwarnings; filterwarnings("ignore")

sc.settings.verbosity = 3

def anndata_copy_attributes(src_anndata, dtn_anndata):
    dtn_anndata.obs = src_anndata.obs
    dtn_anndata.uns = src_anndata.uns
    dtn_anndata.obsm = src_anndata.obsm
    dtn_anndata.obsp = src_anndata.obsp


## collect excel files from all positions
# path = './expression_matrices/210828/220117_analyzed'
path = './expression_matrices/220116'
# path = './expression_matrices/211229'
filepath = os.path.join(path, '*.h5ad')    
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
    adata = sc.concat(adatas, join='inner')
    adata.obs_names_make_unique()

print(adata)
print(adata.var_names)
# print(adata.obs['position'])

print(f'>>> total cell number: {adata.n_obs}')

# %%

adata_virus = adata[:, adata.var_names.isin(virus_list)]
adata_endo = adata[:, ~adata.var_names.isin(virus_list)]

# data quality check
sc.pp.calculate_qc_metrics(adata_endo, percent_top=None, inplace=True, log1p=False)
sc.pl.violin(adata_endo, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)
sc.pl.scatter(adata_endo, x='total_counts', y='n_genes_by_counts')

# %%
print(adata_endo)
# %%
poor_expression_threshold = np.percentile(adata_endo.obs['total_counts'], percentile_filter)
if (poor_expression_threshold == 0) | (np.isnan(poor_expression_threshold)):
    poor_expression_threshold = 1
print(f'poor_expression_threshold: {poor_expression_threshold}')

sns.histplot(adata_endo.obs['total_counts'])

# %%
cell_subset, _ = sc.pp.filter_cells(adata_endo, min_counts=poor_expression_threshold, inplace=False)
cell_subset2, _ = sc.pp.filter_cells(adata_endo, max_counts=2000, inplace=False)
cell_subset = cell_subset & cell_subset2
gene_subset, number_per_gene = sc.pp.filter_genes(adata_endo, min_counts=1, inplace=False)

adata_endo = adata_endo[cell_subset, gene_subset]
adata_virus = adata_virus[cell_subset, :]
print(f'>>> total cells passed the filter: {adata_endo.n_obs}')
sns.histplot(adata_endo.obs['total_counts'])

adata_endo.uns['cell_subset'] = cell_subset

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
gene = 'calb1'
idx = adata_endo.var_names.get_loc(gene)
nonzeros = adata_endo.raw.X[:,idx] > 0
# nonzeros = adata_endo.raw.X[:,idx]>-1
sns.histplot(adata_endo.raw.X[nonzeros,idx], ax=axs[0,0])
sns.histplot(adata_endo_norm.X[nonzeros,idx], ax=axs[0,1])
sns.histplot(adata_endo_log.X[nonzeros,idx], ax=axs[1,0])
sns.histplot(adata_endo_scale.X[nonzeros,idx], ax=axs[1,1])

plt.show()


# %%
sc.pp.normalize_total(adata_endo)
# sc.pp.log1p(adata_endo)
adata_endo.raw = adata_endo
sc.pp.scale(adata_endo)

# %%
for idx in range(adata_endo.n_vars):
    sns.scatterplot(x=adata_endo.obs['total_counts'], y=adata_endo.raw.X[:,idx], y_jitter=0.1)

# %%
for idx in range(adata_endo.n_vars):
    sns.scatterplot(x=adata_endo.obs['total_counts'], y=adata_endo.X[:, idx], y_jitter=0.1)

# %%
sc.pp.pca(adata_endo)
sc.pp.neighbors(adata_endo, use_rep='X_pca', n_neighbors=n_neighbors)
sc.tl.umap(adata_endo)
sc.tl.tsne(adata_endo, use_rep='X_pca')
sc.tl.leiden(adata_endo, resolution=leiden_resolution)
sc.tl.dendrogram(adata_endo, groupby='leiden', use_rep='X_pca', linkage_method='complete')


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

# cluster_reorder = (2, 5, )


# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
print(adata_virus)


# %%

sc.pl.umap(
    adata_virus, 
    color=['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.B8', 'PHP.V1', 'leiden'], 
    s=50, 
    frameon=True, 
    ncols=3, 
    vmax='p99'
)


sc.pl.tsne(
    adata_virus, 
    color=['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.B8', 'PHP.V1', 'leiden'], 
    s=50, 
    frameon=True, 
    ncols=3, 
    vmax='p99'
)

sc.pl.tsne(
    adata_endo,
    color=['slc17a7', 'gad1', 'leiden'],
    # s=50,
    frameon=True,
    vmax='p99'
)

clustering_method='leiden'

# %%
with rc_context({'figure.figsize': (15, 15)}):
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
        vmin=vmin
    )



# %%
# position vs leiden cluster (position-file vs cluster heatmap)
positions = np.asarray(adata_endo.obs['position'])
positions_int = []
for position in positions:
    positions_int.append(int(position[-2:]))
positions = positions_int
labels = np.asarray(adata_endo.obs['leiden']).astype(np.int)
n_pos = np.amax(positions)+1
n_labels = np.amax(labels)+1

position_by_cluster = np.zeros((n_pos, n_labels), dtype=np.int)

for pos, label in zip(positions, labels):
    position_by_cluster[pos, label] = position_by_cluster[pos,label] + 1

# active_position = np.sum(position_by_cluster, axis=1) > 0
active_position = [
    0, 11, 12, 23,
    1, 10, 13, 22,
    2, 9, 14, 21,
    3, 8, 15, 20,
    4, 7, 16, 19,
    5, 6, 17, 18,
]
position_by_cluster = position_by_cluster[np.array(active_position),:]
# active_position = np.argwhere(active_position).flatten()
# print(active_position) 
s = sns.heatmap(data=position_by_cluster, cmap='viridis', yticklabels=active_position)
s.set(xlabel='leiden', ylabel='position')

# %%
total_position_count_per_cluster = np.sum(position_by_cluster, axis=0)

position_by_cluster_norm = position_by_cluster.T/total_position_count_per_cluster[:,None]
position_by_cluster_norm = position_by_cluster_norm.T

s = sns.heatmap(
    data=position_by_cluster_norm,
    cmap='viridis',
    yticklabels=active_position
)
s.set(xlabel='leiden', ylabel='position')

# %%
s = sns.clustermap(
    data=position_by_cluster_norm,
    cmap='viridis',
    yticklabels=active_position,
    row_cluster=False,
    method='ward',
    # vmax=0.5
)
# s.set(xlabel='leiden', ylabel='position')
# s.fig.savefig('./figures/position_cluster.svg', dpi=600)

leiden_idx_reordered_int = s.dendrogram_col.reordered_ind

leiden_idx_reordered_int = [
    9, 10, 6, 0, 18, 17, 15, 13, 5, 12, 16, 2, 3, 14, 7, 1, 4, 8, 11
]
gene_list_ordered = [
    'PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8',
    'slc30a3', 'slc17a7', 'foxp1', 'foxp2', 
    'cux2', 'calb1', 'lamp5',
    'rorb', 'necab1', 'hsd11b1', 'lgi2', 'pcp4', 'crym',
    'rprm', 'tpbg', 'trh', 'ctgf', 
    'calb2', 'crh',
    'ppp1r17',
    'sulf2', 'th', 'chrna6', 'krt73',
    'gad1', 'gad2', 'pvalb', 'reln', 'sst', 'vip', 'sncg', 'tac1'
]



leiden_idx_reordered = [str(idx) for idx in leiden_idx_reordered_int]
adata_endo.obs['leiden_new'] = adata_endo.obs['leiden'].cat.reorder_categories(leiden_idx_reordered)



# %%
# cluster_labels = np.asarray(adata_endo.obs['leiden']).astype(np.uint)
# n_clusters = np.unique(cluster_labels).size
# n_variants = adata_virus.n_vars

# slc17a7_totalcounts = np.zeros((n_clusters,))

# for i, cluster in enumerate(np.unique(cluster_labels)):
#     cluster = int(cluster)
#     cell_inds = np.argwhere(cluster_labels==cluster).ravel()
#     slc17a7_totalcounts[i] = np.mean(adata_endo[cell_inds, 'gad1'].X)

# leiden_idx_reordered_int = list(np.argsort(slc17a7_totalcounts)[::-1])
# leiden_idx_reordered = [str(idx) for idx in leiden_idx_reordered_int]
# adata_endo.obs['leiden_new'] = adata_endo.obs['leiden'].cat.reorder_categories(leiden_idx_reordered)

with rc_context({'figure.figsize': (15, 15)}):
    sc.pl.heatmap(
        adata_endo, 
        # adata_endo.var_names, 
        gene_list_ordered[6:],
        groupby='leiden_new', 
        cmap=cmap_z, 
        # dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=-1
    )
    sc.pl.matrixplot(
        adata_endo, 
        # adata_endo.var_names, 
        gene_list_ordered[6:],
        groupby='leiden_new',
        cmap=cmap_z, 
        # dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin
    )

img = sc.pl.MatrixPlot(
    adata_endo,
    gene_list_ordered[6:],
    groupby='leiden_new',
    vmax=vmax,
    vmin=vmin,
    use_raw=False
).style(cmap_z).swap_axes()
img.savefig('./figures/layer_endo.svg', format='svg')


# %%
cluster_labels = np.asarray(adata_endo.obs['leiden']).astype(np.uint)
n_clusters = np.unique(cluster_labels).size
n_variants = adata_virus.n_vars

transduction_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
transcription_rate = np.zeros((n_variants, n_clusters), dtype=np.uint)
print(np.unique(cluster_labels))
n_cells_in_cluster = np.zeros((n_clusters,))
for i, virus in enumerate(virus_list):
    for j, cluster in enumerate(np.unique(cluster_labels)):
        cluster = int(cluster)
        cell_inds = np.argwhere(cluster_labels==cluster).ravel()
        n_cells_in_cluster[j] = cell_inds.size
        transduction_rate[i,j] = np.count_nonzero(adata_virus[cell_inds, virus].X)*100. / cell_inds.size
        nonzeros = (adata_virus[cell_inds, virus].X != 0).ravel()
        tcpr = np.mean(adata_virus[cell_inds, virus].X[nonzeros].ravel())
        if ~np.isnan(tcpr):
            transcription_rate[i,j] = tcpr

print(n_cells_in_cluster)
transduction_rate = transduction_rate[:, np.array(leiden_idx_reordered_int)]
transcription_rate = transcription_rate[:, np.array(leiden_idx_reordered_int)]

# %%
from matplotlib.pyplot import cm
from matplotlib.colors import LinearSegmentedColormap
blues_zero_white = LinearSegmentedColormap.from_list('', ['white', *cm.Blues(np.arange(255))])
reds_zero_white = LinearSegmentedColormap.from_list('', ['white', *cm.Reds(np.arange(255))])
fig, axs = plt.subplots(2, 1, figsize=(18,9))
sns.heatmap(
    data=transduction_rate, 
    cmap=blues_zero_white, 
    # xticklabels=np.unique(cluster_labels),
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list, 
    ax=axs[0],
    linewidths=.5,
    vmin=0
)
sns.heatmap(
    data=transcription_rate, 
    cmap=reds_zero_white, 
    # xticklabels=np.unique(cluster_labels),
    xticklabels=leiden_idx_reordered,
    yticklabels=virus_list, 
    ax=axs[1],
    linewidths=.5,
    vmin=0
)
fig.savefig('./figures/layer_virus.svg', format='svg')

# %%
# adata_endo.write_h5ad(os.path.join(path, 'endo.h5ad'))
# adata_virus.write_h5ad(os.path.join(path, 'virus.h5ad'))
# # %%
# adata_endo_saved = sc.read_h5ad(os.path.join(path, 'endo.h5ad'))
# print(adata_endo_saved)
