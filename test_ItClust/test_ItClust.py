import ItClust as ic
import scanpy.api as sc
import anndata
import os
from numpy.random import seed
from tensorflow import set_random_seed
import pandas as pd
import numpy as np
import warnings
#os.environ["CUDA_VISIBLE_DEVICES"]="1"
warnings.filterwarnings("ignore")
#import sys
#!{sys.executable} -m pip install 'scanpy==1.4.4.post1'
#Set seeds
#RANDOM_SEED=0
#seed(RANDOM_SEED)
#np.random.seed(RANDOM_SEED)
#set_random_seed(RANDOM_SEED) # on GPU may be some other default


from preprocess import load_PBMC_data

def ItClust_process(data_dir, result_dir):
    ctrl_adata = load_PBMC_data.load_PBMC_batch2_ctrl_data(data_dir)
    stim_adata = load_PBMC_data.load_PBMC_batch2_stim_data(data_dir)
    # use the same control data as target
    #stim_adata = load_PBMC_data.load_PBMC_batch2_ctrl_data(data_dir)
 
    ## load meta data information
    PBMC_df = load_PBMC_data.load_PBMC_batch2_df(os.path.dirname(data_dir))

    ## remain singlet cells
    intersected = list(set(ctrl_adata.obs['barcode']) & set(PBMC_df['barcode'])) 
    sub_ctrl_adata = ctrl_adata[ctrl_adata.obs.barcode.isin(intersected), :]

    intersected = list(set(stim_adata.obs['barcode']) & set(PBMC_df['barcode'])) 
    sub_stim_adata = stim_adata[stim_adata.obs.barcode.isin(intersected), :]

    # add condition
    sub_ctrl_adata.obs['condition'] = 'ctrl'
    sub_stim_adata.obs['condition'] = 'stim'

    # add individual info
    sub_ctrl_adata.obs['ind'] = [PBMC_df[PBMC_df['barcode'] == barcode].iloc[0]['ind'] for barcode in sub_ctrl_adata.obs['barcode']]
    sub_stim_adata.obs['ind'] = [PBMC_df[PBMC_df['barcode'] == barcode].iloc[0]['ind'] for barcode in sub_stim_adata.obs['barcode']]
    sub_ctrl_adata.obs['celltype'] = [PBMC_df[PBMC_df['barcode'] == barcode].iloc[0]['cell'] for barcode in sub_ctrl_adata.obs['barcode']]
    #sub_stim_adata.obs['celltype'] = [PBMC_df[PBMC_df['barcode'] == barcode].iloc[0]['cell'] for barcode in sub_stim_adata.obs['barcode']]

    clf=ic.transfer_learning_clf()
    clf.fit(sub_ctrl_adata, sub_stim_adata)

    pred, prob, cell_type_pred=clf.predict(save_dir=result_dir)
    
    # run t-SNE
    tsne_df = pd.DataFrame(clf.tSNE())
    tsne_df.index = clf.adata_test.obs.index
    tsne_df.to_csv(result_dir+os.sep+"ItClust_tsne.csv")

    return clf

