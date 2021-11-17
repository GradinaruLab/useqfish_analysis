# %%
import scanpy as sc
import numpy as np
import pandas as pd

from glob import glob
import os

from params import *
from matplotlib.pyplot import rc_context

from warnings import filterwarnings; filterwarnings("ignore")

sc.settings.verbosity = 3

def anndata_copy_attributes(src_anndata, dtn_anndata):
    dtn_anndata.obs = src_anndata.obs
    dtn_anndata.uns = src_anndata.uns
    dtn_anndata.obsm = src_anndata.obsm
    dtn_anndata.obsp = src_anndata.obsp

## collect excel files from all positions
# path = './expression_matrices'
path = './expression_matrices/210828'
filepath = os.path.join(path, '*.xlsx')    
filenames = sorted(glob(filepath), key=os.path.basename)
print(filenames)
# filenames = filenames[:3]

resol = 1

datas = [pd.read_excel(filename, index_col=0) for filename in filenames]
for data in datas:
    data.drop(0, inplace=True)
datas = pd.concat(datas, axis=0, ignore_index=True)
gene_dict = {i:gene for i, gene in enumerate(gene_list)}
datas.rename(gene_dict, axis='columns', inplace=True)
datas.drop(columns=[None], inplace=True)
datas = datas[gene_list_ordered]
# print(datas.info())

## convert the concatenated dataframe to anndata
adata = sc.AnnData(datas)
adata_endo = sc.AnnData(datas.iloc[:,n_variants:])
adata_virus = sc.AnnData(datas.iloc[:,:n_variants])
# adata = adata_endo

print(f'>>> total cell number: {adata.n_obs}')


# %%
## preprocessing
cell_subset, number_per_cell = sc.pp.filter_cells(adata_endo, min_genes=1, inplace=False)
gene_subset, number_per_gene = sc.pp.filter_genes(adata_endo, min_counts=1, inplace=False)

adata_endo = adata_endo[cell_subset, gene_subset]
adata_virus = adata_virus[cell_subset, :]
adata = adata[cell_subset, np.concatenate([np.ones((6,), dtype=bool),gene_subset])]

print(f'>>> total cells passed the filter: {adata_endo.n_obs}')

# adata[:, adata.var_names.str.match('genenames')]

# %%
# sc.pl.highest_expr_genes(adata, n_top=10, )

# %%
sc.pp.normalize_total(adata_endo)
sc.pp.log1p(adata_endo)
sc.pp.neighbors(adata_endo)
sc.tl.umap(adata_endo)
sc.tl.tsne(adata_endo)

# %%
sc.tl.leiden(adata_endo, resolution=resol)
sc.tl.dendrogram(adata_endo, groupby='leiden')



# %%
print(adata_endo)
anndata_copy_attributes(adata_endo, adata_virus)
anndata_copy_attributes(adata_endo, adata)
print(adata_virus)
sc.pp.log1p(adata_virus)

# %%
# sc.pl.umap(adata)
# sc.pl.embedding(adata_endo, basis="X_umap", color="leiden")


# %%
sc.tl.rank_genes_groups(adata_endo, groupby='leiden')

# %%


with rc_context({'figure.figsize': (3, 3)}):
    sc.pl.umap(
        adata_virus, 
        color=['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.B8', 'PHP.V1', 'leiden'], 
        s=50, 
        frameon=True, 
        ncols=3, 
        vmax='p99'
    )

with rc_context({'figure.figsize': (3, 3)}):
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
sc.pl.heatmap(
    adata_endo, 
    adata_endo.var_names, 
    groupby='leiden', 
    cmap='viridis', 
    # dendrogram=True,
    swap_axes=True
)
sc.pl.heatmap(
    adata_virus, 
    adata_virus.var_names, 
    groupby='leiden', 
    cmap='viridis', 
    # dendrogram=True,
    swap_axes=True
)

# %%
viral_gene_dict = {gene_name: [gene_name] for gene_name in ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.B8', 'PHP.V1']}
endo_gene_dict = {gene_name: [gene_name] for gene_name in gene_list_ordered}
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

sc.pl.matrixplot(
    adata_endo,
    adata_endo.var_names,
    'leiden',
    swap_axes=True

)
sc.pl.matrixplot(
    adata_virus, 
    # viral_gene_dict, 
    adata_virus.var_names,
    'leiden', 
    cmap='Blues', 
    standard_scale='var', 
    colorbar_title='column scaled\nexpression',
    swap_axes=True,
)
sc.pl.matrixplot(
    adata_virus, 
    # viral_gene_dict, 
    adata_virus.var_names,
    'leiden',
    swap_axes=True
)
# sc.pl.dotplot(adata_endo, marker_genes_dict, 'leiden')

# %%
sc.tl.rank_genes_groups(adata_endo, 'leiden')
sc.pl.rank_genes_groups(adata_endo, n_genes=5, sharey=False)


# # %%
# sc.pp.normalize_total(adata_virus)
# sc.pp.log1p(adata_virus)
# sc.pp.neighbors(adata_virus)
# sc.tl.umap(adata_virus)

# # %%
# sc.tl.leiden(adata_virus, resolution=1)
# # sc.tl.dendrogram(adata_virus, groupby='leiden')



# # %%
# # print(adata_endo)
# anndata_copy_attributes(adata_virus, adata_endo)
# anndata_copy_attributes(adata_virus, adata)
# # print(adata_virus)
# sc.pp.log1p(adata_endo)

# # %%
# # sc.pl.umap(adata)
# # sc.pl.embedding(adata_endo, basis="X_umap", color="leiden")


# # # %%
# # sc.tl.rank_genes_groups(adata_endo, groupby='cluster')
# # sc.pl.rank_genes_groups(adata_endo, sharey=False)


# # %%
# sc.pl.heatmap(
#     adata_virus, 
#     adata_virus.var_names, 
#     groupby='leiden', 
#     cmap='viridis', 
#     swap_axes=True
# )
# sc.pl.heatmap(
#     adata_endo, 
#     adata_endo.var_names, 
#     groupby='leiden', 
#     cmap='viridis', 
#     swap_axes=True
# )

# with rc_context({'figure.figsize': (3, 3)}):
#     sc.pl.umap(
#         adata_virus, 
#         color=['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.B8', 'PHP.V1', 'leiden'], 
#         s=50, 
#         frameon=True, 
#         ncols=2, 
#         vmax='p99'
#     )
# # %%

# %%
