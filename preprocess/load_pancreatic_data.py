import os
import anndata
import numpy as np
import pandas as pd

# --------------
# load pancreatic counts and labels from Ziyi's data
# --------------
def load_panc_data(data_dir, source="seg"):
    '''Loading pancreatic islets into anndata
    #@source: seg/muraro/xin
    '''
    ## load data as dataframe
    if "seg" == source:
        df = pd.read_csv(data_dir+os.sep+"Pancreas6/seg/Seg_counts.csv", index_col=0)
        obs = pd.read_csv(data_dir+os.sep+"Pancreas6/seg/Seg_label.csv", index_col=0)
    elif "muraro" == source:
        df = pd.read_csv(data_dir+os.sep+"Pancreas6/muraro/Muraro_counts.csv", index_col=0)
        obs = pd.read_csv(data_dir+os.sep+"Pancreas6/muraro/Muraro_label.csv", index_col=0)
    elif "xin" == source:
        df = pd.read_csv(data_dir+os.sep+"Pancreas6/xin/Xin_counts.csv", index_col=0)
        obs = pd.read_csv(data_dir+os.sep+"Pancreas6/xin/Xin_label.csv", index_col=0)
    ## build obs/var dataframe
    obs.columns = ["cell.type"]
    obs.index = df.columns
    var = pd.DataFrame(data=df.index, index=df.index)
    var.columns = ['gene_symbols']

    adata = anndata.AnnData(X=df.T, obs=obs, var=var)
    adata.obs_names_make_unique(join="-")
    adata.var_names_make_unique(join="-")
    return adata

def load_seg_data(data_dir, condition="Healthy"):
    '''Load seg dataset and split into two conditions
    #@condition: Healthy/T2D
    '''
    adata = load_panc_data(data_dir, source="seg")

    ## add healthy/T2D condition to 
    adata.obs['condition'] = ["T2D" if "T2D" in cell else "Healthy" for cell in adata.obs.index]
    condition_cells = adata.obs[adata.obs['condition'] == condition].index
    return adata[condition_cells]

def load_xin_data(data_dir, condition="Healthy"):
    '''Load xin dataset and split into two conditions
    #@condition: Healthy/T2D
    '''
    adata = load_panc_data(data_dir, source="xin")

    #annots = 
    return adata
