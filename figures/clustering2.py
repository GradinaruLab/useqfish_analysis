
# %% 

# parameters
n_neighbors = 20
vmax = 3
vmin = -1
cmap_z = 'coolwarm'
cmap_raw = 'viridis'
poor_cluster_percentile_threshold = 50
poor_cluster_zstd_threshold = 10
visualization = True
leiden_resolution = .5


# %%

# data loading
import scanpy as sc
import scanpy.external as scex
import seaborn as sns

import numpy as np
from scipy.stats import mannwhitneyu

from glob import glob
import os

from params import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

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
# path = './expression_matrices/210828/220122_analyzed'
path = './expression_matrices/major'
# path = './expression_matrices/h5ad_old'
# path = './expression_matrices/210828'
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

adata.obs['cell_id'] = pd.Categorical(values=range(adata.n_obs))


print(adata)
print(adata.var_names)
# print(adata.obs['position'])
# adata.write(os.path.join(path, 'expression_matrix.h5ad'))

print(f'>>> total cell number: {adata.n_obs}')

# %%
# endogenous and viral gene separation
adata_virus = adata[:, adata.var_names.isin(virus_list)]
adata_endo = adata[:, ~adata.var_names.isin(virus_list)]

# data quality check
sc.pp.calculate_qc_metrics(adata_endo, percent_top=None, inplace=True, log1p=False)
sc.pl.violin(adata_endo, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)
sc.pl.scatter(adata_endo, x='total_counts', y='n_genes_by_counts')

# %%
# data filtering
poor_expression_threshold = np.percentile(adata_endo.obs['total_counts'], 20)
if poor_expression_threshold == 0:
    poor_expression_threshold = 1
print(f'poor_expression_threshold: {poor_expression_threshold}')

sns.histplot(adata_endo.obs['total_counts'])

# %%
cell_subset, _ = sc.pp.filter_cells(adata_endo, min_counts=poor_expression_threshold, inplace=False)
cell_subset2, _ = sc.pp.filter_cells(adata_endo, max_counts=1e3, inplace=False)
cell_subset = cell_subset & cell_subset2
gene_subset, number_per_gene = sc.pp.filter_genes(adata_endo, min_counts=3, inplace=False)

adata_endo = adata_endo[cell_subset, gene_subset]
adata_virus = adata_virus[cell_subset, :]
print(f'>>> total cells passed the filter: {adata_endo.n_obs}')

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
nonzeros = adata_endo.raw.X[:,idx]>-1
sns.scatterplot(x=adata_endo.obs['total_counts'], y=adata_endo.raw.X[:,idx], y_jitter=0.1, ax=axs[0,0])
sns.histplot(adata_endo.raw.X[nonzeros,idx], ax=axs[0,1])
sns.histplot(adata_endo_norm.X[nonzeros,idx], ax=axs[1,0])
sns.histplot(adata_endo_log.X[nonzeros,idx], ax=axs[1,1])
sns.histplot(adata_endo_scale.X[nonzeros,idx], ax=axs[2,0])
sns.scatterplot(x=adata_endo.obs['total_counts'], y=adata_endo_scale.X[:,idx], y_jitter=0.1, ax=axs[2,1])

plt.show()

# %%
sc.pp.normalize_total(adata_endo)
# sc.pp.log1p(adata_endo)
sc.pp.scale(adata_endo)

# %%
for idx in range(adata_endo.n_vars):
# idx = 11
    sns.scatterplot(x=adata_endo.obs['total_counts'], y=adata_endo.raw.X[:,idx], y_jitter=0.1)

# %%
for idx in range(adata_endo.n_vars):
    sns.scatterplot(x=adata_endo.obs['total_counts'], y=adata_endo.X[:, idx], y_jitter=0.1)

# %%
# data clustering
sc.pp.pca(adata_endo)
sc.pl.pca_overview(adata_endo)
sc.pp.neighbors(adata_endo, use_rep='X_pca', n_neighbors=n_neighbors)
sc.tl.umap(adata_endo)
sc.tl.tsne(adata_endo, use_rep='X_pca')
sc.tl.leiden(adata_endo, resolution=leiden_resolution)
sc.tl.dendrogram(adata_endo, groupby='leiden', use_rep='X_pca', linkage_method='ward')

# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
print(adata_virus)

#  %%
# clustering visualization
sc.pl.umap(
    adata_virus, 
    color=['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.B8', 'PHP.V1', 'leiden'], 
    frameon=True, 
    ncols=3, 
    vmax='p99'
)
sc.pl.tsne(
    adata_virus, 
    color=['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.B8', 'PHP.V1', 'leiden'], 
    frameon=True, 
    ncols=3, 
    vmax='p99'
)


clustering_method='leiden'
   
# %%
# heatmaps
with rc_context({'figure.figsize': (15, 15)}):
    sc.pl.heatmap(
        adata_endo, 
        ['slc17a7', 'gad1', 'pvalb', 'sst', 'vip', 'gja1', 'mbp', 'msr1', 'hexb', 'cldn5', 'acta2'], 
        groupby=clustering_method, 
        cmap=cmap_z, 
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin
        # log=True
    )
    sc.pl.matrixplot(
        adata_endo,
        ['slc17a7', 'gad1', 'pvalb', 'sst', 'vip', 'gja1', 'mbp', 'msr1', 'hexb', 'cldn5', 'acta2'],
        groupby=clustering_method,
        swap_axes=True,
        # dendrogram=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin,
        cmap=cmap_z
    )

# %%
# subclustering
adata_sub = adata_endo[adata_endo.obs['leiden']=='1', :]
print(adata_sub)
sc.pp.pca(adata_sub)
sc.pp.neighbors(adata_sub, use_rep='X_pca')
sc.tl.leiden(adata_sub, key_added='leiden_sub', resolution=leiden_resolution)
sc.tl.dendrogram(adata_sub, groupby='leiden_sub', linkage_method='ward')

sc.pl.heatmap(
    adata_sub,
    ['slc17a7', 'gad1', 'pvalb', 'sst', 'vip', 'gja1', 'mbp', 'msr1', 'hexb', 'cldn5', 'acta2'], 
    groupby='leiden_sub',
    cmap=cmap_z,
    dendrogram=True,
    swap_axes=True,
    use_raw=False,
    vmax=vmax,
    vmin=vmin
)

# %%
print(adata_sub)

cluster_labels = np.asarray(adata_endo.obs['leiden']).astype(np.uint)
n_clusters = np.unique(cluster_labels).size
new_cluster_labels = []

# # print(adata_endo[5].obs['leiden'].index)
# print(adata_endo.obs['leiden'].loc['17'])
# print(adata_sub[0].obs['leiden_sub'])

# print(adata_endo[0].obs['leiden'])
# adata_endo.obs['leiden_new'].loc[adata_endo[0].obs['leiden'].index] = adata_endo[0].obs['leiden']
# print(adata_endo[0].obs['leiden_new'])
# print(adata_sub.obs['leiden_sub'])

for cell in range(adata_endo.n_obs):
    # print(adata_endo[cell].obs['leiden'])    
    if adata_endo[cell].obs['leiden'][0] != '1':
        new_cluster_labels.append(adata_endo[cell].obs['leiden'][0])
    else:
        cell_id = adata_endo[cell].obs['leiden'].index
        subcluster_label = adata_sub.obs['leiden_sub'].loc[cell_id][0]
        new_cluster_labels.append(str(n_clusters+1+int(subcluster_label)))
        
        # subind = (adata_sub.obs['cell_id'][0] == cell)
        # break
        # print(f'subind: {subind}')
        # leiden_new.append(n_clusters+1+adata_sub[subind].obs['leiden'][0])
    # print(f'cell: {cell}, {new_cluster_labels}')


adata_endo.obs['leiden_sub'] = pd.Categorical(new_cluster_labels)
print(adata_endo.obs['leiden_sub'])
cluster_list = [str(idx) for idx in range(adata_endo.obs['leiden_sub'].nunique())]
print(cluster_list)

# %%
adata_endo.obs['leiden_sub'] = adata_endo.obs['leiden_sub'].cat.rename_categories(cluster_list)

# %%
# clustering adjustment using dendrogram 
# - compare close clusters with two sided Mann-Whitney U test
# and if not significant, merge two clusters
sc.tl.dendrogram(adata_endo, groupby='leiden_sub', use_rep='X_pca')
with rc_context({'figure.figsize': (18, 15)}):
    s = sc.pl.heatmap(
        adata_endo,
        ['slc17a7', 'gad1', 'pvalb', 'sst', 'vip', 'gja1', 'mbp', 'msr1', 'hexb', 'cldn5', 'acta2'], 
        groupby='leiden_sub',
        cmap=cmap_z,
        dendrogram=True,
        # swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin
    )


# %%
print(adata_endo.uns['dendrogram_leiden_sub']['linkage'])

linkage_matrix = adata_endo.uns['dendrogram_leiden_sub']['linkage']
# print(linkage_matrix.shape)

cluster_labels = np.asarray(adata_endo.obs['leiden_sub']).astype(np.uint)
n_clusters = np.unique(cluster_labels).size
new_labels = {str(i):str(i) for i in range(n_clusters)}
p_threshold = 1e-3

target_n_clusters = 11
for iter in range(linkage_matrix.shape[0] - target_n_clusters+1):
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
    cluster_labels[np.argwhere(cluster_labels==cluster1).ravel()] = n_clusters
    cluster_labels[np.argwhere(cluster_labels==cluster2).ravel()] = n_clusters
    # print(cluster_labels)
    
    # if np.any(res.pvalue > p_threshold):
    # if linkage_matrix[iter,2] < 1:
    # print(list(new_labels.values()).index(str(cluster1)))
    
    cinds = [i for i, cluster in new_labels.items() if cluster==str(cluster1)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    cinds = [i for i, cluster in new_labels.items() if cluster==str(cluster2)]
    for cind in cinds:
        new_labels[cind] = str(n_clusters)
    # new_labels[cinds] = str(n_clusters)
    # new_labels[str(list(new_labels.values()).index(str(cluster2)))] = str(n_clusters)
    # new_labels[str(cluster2)] = str(n_clusters)
    # print(f'clusters merged: {cluster1}+{cluster2} -> {n_clusters}')
    # print(new_labels)
    
    n_clusters = n_clusters + 1

# print(adata_endo.obs['leiden'])
# print(new_labels)
# %%
# print(adata_endo.obs['leiden_sub'])
# adata_endo.obs['leiden_new'] = adata_endo.obs['leiden_sub'].map(new_labels).astype('category')
adata_endo.obs['leiden_new'] = adata_endo.obs['leiden_sub'].map(new_labels).astype('category')
print(adata_endo.obs['leiden_new'])
sc.tl.dendrogram(adata_endo, groupby='leiden_new', use_rep='X_pca')

# %%
# print(adata_endo.obs['leiden_new'])


# %%
# copy attrs acquired from endo to virus
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
print(adata_virus)



# %%
# clustering visualization
clustering_method='leiden_new'
sc.pl.umap(
    adata_virus, 
    color=['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.B8', 'PHP.V1', clustering_method], 
    frameon=True, 
    ncols=3, 
    vmax='p99'
)
sc.pl.tsne(
    adata_virus, 
    color=['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.B8', 'PHP.V1', clustering_method], 
    frameon=True, 
    ncols=3, 
    vmax='p99'
)

# heatmaps
with rc_context({'figure.figsize': (15, 15)}):
    sc.pl.heatmap(
        adata_endo, 
        ['slc17a7', 'gad1', 'pvalb', 'sst', 'vip', 'gja1', 'mbp', 'msr1', 'hexb', 'cldn5', 'acta2'], 
        groupby=clustering_method, 
        cmap=cmap_z, 
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin
    )
    sc.pl.matrixplot(
        adata_endo,
        ['slc17a7', 'gad1', 'pvalb', 'sst', 'vip', 'gja1', 'mbp', 'msr1', 'hexb', 'cldn5', 'acta2'],
        groupby=clustering_method,
        swap_axes=True,
        # dendrogram=True,
        use_raw=False,
        vmax=vmax,
        vmin=vmin,
        cmap=cmap_z,
    )


# %%
# clustering_method='leiden'
cluster_labels = np.asarray(adata_endo.obs[clustering_method]).astype(np.uint)
print(cluster_labels)
n_clusters = np.unique(cluster_labels).size
print(n_clusters)
n_variants = adata_virus.n_vars

transduction_rate = np.zeros((n_variants, n_clusters))
transcription_rate = np.zeros((n_variants, n_clusters))

for i, virus in enumerate(virus_list):
    for j, cluster in enumerate(np.unique(cluster_labels)):
        cluster = int(cluster)
        cell_inds = np.argwhere(cluster_labels==cluster).ravel()
        transduction_rate[i,j] = np.count_nonzero(adata_virus[cell_inds, virus].X)*100. / cell_inds.size
        nonzeros = (adata_virus[cell_inds, virus].X != 0)
        transcription_rate[i,j] = np.mean(adata_virus[cell_inds, virus].X[nonzeros])


transcription_rate_scaled = (transcription_rate - transcription_rate.min(axis=1)[:,None]) / transcription_rate.max(axis=1)[:,None]

leiden_sorted_order = [
    '19', '20', '5', '17', '14', '18', '15', '12', '4', '3', '2'
]

leiden_new_order = list(np.unique(np.asarray(adata_endo.obs[clustering_method]).astype(np.uint)))
print(leiden_new_order)
index = [leiden_new_order.index(int(sortedind)) for sortedind in leiden_sorted_order]

print(index)

# %%
fig, axs = plt.subplots(2, 1, figsize=(18,9))
sns.heatmap(
    data=transduction_rate[:, index], 
    cmap='viridis', 
    xticklabels=leiden_sorted_order,
    yticklabels=virus_list, 
    ax=axs[0],
    # annot=True
)
sns.heatmap(
    data=transcription_rate[:, index], 
    cmap=reds_zero_white, 
    xticklabels=leiden_sorted_order,
    yticklabels=virus_list, 
    ax=axs[1],
    annot=True
)

# %%
adata_virus2 = adata_virus
# adata_virus2.X = adata_virus2.X/adata_virus2.X.sum(axis=0)[None,:]*adata_virus2.n_vars
# adata_virus2.X = adata_virus2.X/adata_virus2.X.sum(axis=1)[:,None]*adata_virus2.n_obs

# %%
sc.pp.normalize_total(adata_virus)
print(adata_virus.to_df().head())
sc.pp.log1p(adata_virus)
print(adata_virus.to_df().head())
sc.pp.scale(adata_virus)
print(adata_virus.to_df().head())


img = sc.pl.MatrixPlot(
    adata_virus,
    virus_list,
    'leiden_new',
    categories_order=leiden_sorted_order,
    use_raw=False,
    # vmax=.5,
    # vmin=-.5
).style(cmap=cmap_z).swap_axes()
img.show()


# %%
## image save
from matplotlib.pyplot import cm
from matplotlib.colors import LinearSegmentedColormap
blues_zero_white = LinearSegmentedColormap.from_list('', ['white', *cm.Purples(np.arange(255))])
reds_zero_white = LinearSegmentedColormap.from_list('', ['white', *cm.Greens(np.arange(255))])
# cb_orange = LinearSegmentedColormap.from_list('', ['white', '#137837'])
# cb_teal = LinearSegmentedColormap.from_list('', ['white', '#762A83'])

img = sc.pl.MatrixPlot(
    adata_endo,
    ['slc17a7', 'gad1', 'pvalb', 'sst', 'vip', 'gja1', 'mbp', 'hexb', 'msr1', 'cldn5', 'acta2'], 
    'leiden_new',
    categories_order=leiden_sorted_order,
    use_raw=False,
    vmax=vmax,
    vmin=vmin,
    # standard_scale='var'
).style(cmap='summer').swap_axes()
img.savefig('./figures/major_endo.svg', format='svg')

img = sc.pl.MatrixPlot(
    adata_virus,
    virus_list,
    'leiden_new',
    categories_order=leiden_sorted_order,
    use_raw=False,
    vmax=.5,
    vmin=-.5
).style(cmap=cmap_z).swap_axes()
img.savefig('./figures/major_bias.svg', format='svg')

fig, axs = plt.subplots(2, 1, figsize=(18,9))
sns.heatmap(
    data=transduction_rate[:,index],
    cmap=blues_zero_white,
    # cmap=cb_orange,
    xticklabels=leiden_sorted_order,
    yticklabels=virus_list,
    ax=axs[0],
    linewidths=.5,
    vmin=0
)
sns.heatmap(
    data=transcription_rate[:,index],
    cmap=reds_zero_white,
    xticklabels=leiden_sorted_order,
    yticklabels=virus_list,
    ax=axs[1],
    linewidths=.5,
    vmin=0
)
# sns.heatmap(
#     data=transcription_rate_scaled[:,index],
#     cmap='viridis',
#     xticklabels=leiden_sorted_order,
#     yticklabels=virus_list,
#     ax=axs[2],
#     linewidths=.5,
# )
fig.savefig('./figures/major_virus.svg', format='svg')# %%


# %%
