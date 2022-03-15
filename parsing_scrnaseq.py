# %%
import scanpy as sc
import scanpy.external as scex
import numpy as np
import pandas as pd

from glob import glob
import os

from params import *
from matplotlib.pyplot import rc_context

from converting_anndata import *

from warnings import filterwarnings

filterwarnings("ignore")

sc.settings.verbosity = 3


def anndata_copy_attributes(src_anndata, dtn_anndata):
    dtn_anndata.obs = src_anndata.obs
    dtn_anndata.uns = src_anndata.uns
    dtn_anndata.obsm = src_anndata.obsm
    dtn_anndata.obsp = src_anndata.obsp


## collect excel files from all positions
# path = './expression_matrices/211024'
# path = './expression_matrices/211229'
# path = './expression_matrices/220116'
path = "./expression_matrices/scrnaseq"
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
    adata = sc.concat(adatas, join="outer")
    adata.obs_names_make_unique()

print(adata)
print(adata.var_names)

print(f">>> total cell number: {adata.n_obs}")

# %%
ref_gene_list = ["ENSMUSG00000002930"]
adata_ref_gene = adata[:, adata.var_names.isin(ref_gene_list)]
data_ref_gene = adata_ref_gene.X.toarray()
data_ref_gene = data_ref_gene[~np.isnan(data_ref_gene)]
print(f"nonzero_mean_expression: {np.mean(data_ref_gene[np.nonzero(data_ref_gene)])}")
print(f"mean_expression: {np.mean(data_ref_gene)}")


# %%
