import scanpy as sc
import pandas as pd

from glob import glob
import os

from params import *

from warnings import filterwarnings; filterwarnings("ignore")


def xlsx2h5ad(path):
    filepath = os.path.join(path, '*.xlsx')    
    filenames = sorted(glob(filepath), key=os.path.basename)
    print(filenames)

    filepath_adata = os.path.join(path, 'expression_matrix.h5ad')

    datas = [pd.read_excel(filename, index_col=0) for filename in filenames]
    obs_position = []    
    for data, filename in zip(datas, filenames):
        data.drop(0, inplace=True)

        position = filename.split('/')[-1]
        position = position.split('_')[0]
        for i in range(data.shape[0]):
            obs_position.append(position)

    datas = pd.concat(datas, axis=0, ignore_index=True)
    gene_dict = {i:gene for i, gene in enumerate(gene_list)}
    datas.rename(gene_dict, axis='columns', inplace=True)
    # datas.drop(columns=[None], inplace=True)
    datas = datas[gene_list_ordered]
    print(datas.info())

    ## convert the concatenated dataframe to anndata
    adata = sc.AnnData(datas, var=pd.DataFrame(index=gene_list_ordered))
    adata.obs['position'] = obs_position
    
    # adata.X = csr_matrix(adata.X)
    adata.write(filepath_adata)
    print(f'>>>> expression matrices were converted to h5ad format')

