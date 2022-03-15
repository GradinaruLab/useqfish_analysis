# %%
import scanpy as sc
import scanpy.external as scex
import seaborn as sns

import numpy as np

from glob import glob
import os

from params import *
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
# path = './expression_matrices/210828/220117_analyzed'
path = './expression_matrices/h5ad_old'
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
total_counts_per_cell = np.sum(adata_endo.X, axis=1)
poor_expression_threshold = np.percentile(total_counts_per_cell, 50)
if poor_expression_threshold == 0:
    poor_expression_threshold = 1

n_genes_per_cell = np.sum(adata_endo.X.astype(bool), axis=1)
# cell_outliers_threshold = np.percentile(n_genes_per_cell, 70)

# sns.histplot(total_counts_per_cell)
# sns.histplot(n_genes_per_cell)
print(f'poor_expression_threshold: {poor_expression_threshold}')
# print(f'cell_outliers_threshold: {cell_outliers_threshold}')

cell_subset, number_per_cell = sc.pp.filter_cells(adata_endo, min_counts=poor_expression_threshold, inplace=False)
# cell_subset, number_per_cell = sc.pp.filter_cells(adata_endo, min_genes=1, inplace=False)
gene_subset, number_per_gene = sc.pp.filter_genes(adata_endo, min_counts=1, inplace=False)


adata_endo = adata_endo[cell_subset, gene_subset]
adata_virus = adata_virus[cell_subset, :]
# adata = adata[cell_subset, np.concatenate([np.ones((6,), dtype=bool),gene_subset])]

print(f'>>> total cells passed the filter: {adata_endo.n_obs}')

# %%
adata_endo_sim = scex.pp.scrublet_simulate_doublets(adata_endo)
adata_endo = scex.pp.scrublet(
    adata_endo, 
    adata_sim=adata_endo_sim,
    # n_prin_comps=3, 
    n_prin_comps=10,
    copy=True
)
cell_singlet = ~adata_endo.obs['predicted_doublet'].fillna(False).astype(bool)
adata_endo = adata_endo[cell_singlet]
adata_virus = adata_virus[cell_singlet]

# adata_endo = sc.AnnData(datas.iloc[:,n_variants:])
# adata_virus = sc.AnnData(datas.iloc[:,:n_variants])

print(adata_endo)
print(adata_virus)
print(f'>>> total singlet cells: {adata_endo.n_obs}')

# %%
# adata_endo.raw = adata_endo

sc.pp.normalize_total(adata_endo)
# sc.pp.log1p(adata_endo)
adata_endo.raw = adata_endo
print(adata_endo.X.max(), adata_endo.X.min())
# sc.pp.scale(adata_endo)
print(adata_endo.X.max(), adata_endo.X.min())
# sc.pp.highly_variable_genes(adata_endo)
# sc.pp.pca(adata_endo)

# %%
# sc.pp.neighbors(adata_endo, use_rep='X_pca')
sc.pp.neighbors(adata_endo)
sc.tl.umap(adata_endo)
# sc.tl.tsne(adata_endo, use_rep='X_pca')
sc.tl.tsne(adata_endo)

# %%
sc.tl.leiden(adata_endo, resolution=leiden_resolution)
sc.tl.dendrogram(adata_endo, groupby='leiden')


# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
# anndata_copy_attributes(adata_endo, adata)
print(adata_virus)
# sc.pp.log1p(adata_virus)


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
# sc.pl.umap(adata_virus, color='PHP.eB')
# sc.pl.umap(adata_virus, color='CAP-B10')
# sc.pl.umap(adata_virus, color='PHP.N')
# sc.pl.umap(adata_virus, color='PHP.Astro')

# %%
# with rc_context({'figure.figsize': (6, 3)}):
#     for name in adata_virus.var_names:
#         sc.pl.violin(
#             adata_virus, 
#             [name],
#             groupby='leiden',
#         )
    
# %%
with rc_context({'figure.figsize': (15, 15)}):
    sc.pl.heatmap(
        adata_endo, 
        adata_endo.var_names, 
        groupby='leiden', 
        cmap='bwr', 
        dendrogram=True,
        swap_axes=True,
        use_raw=False,
        vmax=2,
        vmin=-2
        # log=True
    )
    sc.pl.heatmap(
        adata_endo,
        adata_endo.var_names, 
        groupby='leiden', 
        cmap='viridis', 
        dendrogram=True,
        swap_axes=True,
    )
    sc.pl.heatmap(
        adata_virus, 
        adata_virus.var_names, 
        groupby='leiden', 
        cmap='viridis', 
        dendrogram=True,
        swap_axes=True,
    )

# cluster_annotation = {
#     '0': 'slc17a7',
#     '1': 'excitatory-like',
#     '2': 'gad1/pvalb',
#     '3': 'sst',
#     '4': 'gja1',
#     '5': 'vip',
#     '6': 'hexb',
#     '7': 'cldn5',
#     '8': 'acta2'
# }
# adata_endo.obs['cell type'] = adata_endo.obs['leiden'].map(cluster_annotation).astype('category')
# adata_virus.obs['cell type'] = adata_virus.obs['leiden'].map(cluster_annotation).astype('category')

# %%
print(adata_endo.var_names)
with rc_context({'figure.figsize': (15, 15)}):
    sc.pl.matrixplot(
        adata_endo,
        # adata_endo.var_names,
        # gene_list_ordered[6:],
        ['slc17a7', 'gad1', 'pvalb', 'sst', 'vip', 'gja1', 'mbp', 'msr1', 'hexb', 'cldn5', 'acta2'],
        'leiden',
        swap_axes=True,
        dendrogram=True,
        use_raw=False,
        vmax=1,
        vmin=-1,
        cmap='bwr'
    )
    sc.pl.matrixplot(
        adata_endo,
        # adata_endo.var_names,
        # gene_list_ordered[6:],
        ['slc17a7', 'gad1', 'pvalb', 'sst', 'vip', 'gja1', 'mbp', 'msr1', 'hexb', 'cldn5', 'acta2'],
        'leiden',
        swap_axes=True,
        dendrogram=True,
        cmap='viridis'
    )
    
    sc.pl.matrixplot(
        adata_virus,
        adata_virus.var_names,
        'leiden',
        cmap='Blues',
        standard_scale='obs',
        colorbar_title='column log scaled\nexpression',
        swap_axes=True,
        dendrogram=True,
    )
    sc.pl.matrixplot(
        adata_virus, 
        adata_virus.var_names,
        'leiden', 
        cmap='Reds', 
        standard_scale='var', 
        colorbar_title='row log scaled\nexpression',
        swap_axes=True,
        dendrogram=True,
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

active_position = np.sum(position_by_cluster, axis=1) > 0
position_by_cluster = position_by_cluster[active_position,:]
active_position = np.argwhere(active_position).flatten()
print(active_position) 
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
cluster_model = AgglomerativeClustering(n_clusters=4).fit(position_by_cluster_norm.T)
position_predicted = cluster_model.labels_
idx = np.argsort(position_predicted)
print(idx)
position_by_cluster_norm_sorted = position_by_cluster_norm[:, idx]
s = sns.heatmap(data=position_by_cluster_norm_sorted, cmap='viridis', xticklabels=idx)

leiden_idx = [str(i) for i in idx]
print(adata_endo.var_names)
with rc_context({'figure.figsize': (15, 15)}):    
    img = sc.pl.MatrixPlot(
        adata_endo,
        # adata_endo.var_names,
        gene_list_ordered[6:],
        'leiden',
        categories_order=leiden_idx,
        use_raw=False,
        vmax=1,
        vmin=-1,
        cmap='bwr'
    ).swap_axes()
    img.show()
    # img.savefig('./figures/cortex_endo.svg', format='svg')

    img = sc.pl.MatrixPlot(
        adata_virus,
        ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8'],
        'leiden',
        categories_order=leiden_idx,
        standard_scale='obs',
    ).style(cmap='Blues').swap_axes()
    img.show()
    # img.savefig('./figures/cortex_virus_obs.svg', format='svg')

    img = sc.pl.MatrixPlot(
        adata_virus,
        ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8'],
        'leiden',
        categories_order=leiden_idx,
        standard_scale='var',
    ).style(cmap='Reds').swap_axes()
    img.show()
    # img.savefig('./figures/cortex_virus_var.svg', format='svg')
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
