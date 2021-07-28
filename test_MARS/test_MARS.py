## should use mars environment

import torch
import os, sys, time
import numpy as np
import pandas as pd
import scanpy.api as sc
from anndata import AnnData
import anndata

from matplotlib import pyplot as plt
import matplotlib as mpl

sys.path.append("../")
sys.path.append('/homelocal/wma36/celltyping_refConstruct/test_MARS/mars')
from args_parser import get_parser
from model.mars import MARS
from model.experiment_dataset import ExperimentDataset
import warnings
warnings.filterwarnings('ignore')

## set seed
torch.manual_seed(0)
import random
random.seed(0)
np.random.seed(0)

## pytorch params
params, unknown = get_parser().parse_known_args()
params.cuda = True
device = 'cuda:0' if torch.cuda.is_available() and params.cuda else 'cpu'
params.device = device

def run_MARS(train_adata, test_adata, result_dir, prefix="MARS", 
        celltype_col="cell.type"):
    '''Run MARS pipeline
    '''

    adata = train_adata.concatenate(test_adata)

    from sklearn import preprocessing
    le = preprocessing.LabelEncoder()
    le.fit(train_adata.obs[celltype_col])
    y_train = np.array(le.transform(train_adata.obs[celltype_col]), dtype=np.int64)
    annotated = ExperimentDataset(train_adata.X, train_adata.obs_names, 
            train_adata.var_names, '0', y_train)
    
    y_test = np.array(le.transform(test_adata.obs[celltype_col]), dtype=np.int64)
    unannnotated = ExperimentDataset(test_adata.X, test_adata.obs_names, 
            test_adata.var_names, '1', y_test)
    
    pretrain_data = ExperimentDataset(test_adata.X, test_adata.obs_names, 
            test_adata.var_names, '1')
    n_clusters = len(np.unique(unannnotated.y))
    
    #params.epochs = 100 ## change to 50/100 and see the difference
    start = time.time()
    mars = MARS(n_clusters, params, [annotated], unannnotated, pretrain_data, hid_dim_1=1000, hid_dim_2=100)
    mars_adata, landmarks, scores = mars.train(evaluation_mode=True, save_all_embeddings=True) # evaluation mode
    end = time.time()
    with open(result_dir+os.sep+prefix+'_metrics.txt', 'w') as f:
        f.write("Acc:%s\n" % str(scores['accuracy']))
        f.write("ARI:%s\n" % str(scores['adj_rand']))
        f.write("macroF1:%s\n" % str(scores['f1_score']))
        f.write("runtime:%s\n" % str(end-start))
    
    ## give names to cell types (Aborted, very inaccurate)
    #keys = le.transform(le.classes_)
    #values = le.classes_
    #cell_type_name_map = dict()
    #for index, key in enumerate(keys):
    #    cell_type_name_map[key] = values[index]
    #interp_names = mars.name_cell_types(mars_adata, landmarks, cell_type_name_map)

    mars_adata.obs = mars_adata.obs.merge(adata.obs, 
            left_index=True, right_index=True, how="left")
    mars_adata = mars_adata[mars_adata.obs['dataset_batch'] == "test"]

    adata_mars = AnnData(mars_adata.obsm['MARS_embedding'])
    adata_mars.obs['MARS_labels'] = pd.Categorical(mars_adata.obs['MARS_labels'])
    adata_mars.obs['ground_truth'] = pd.Categorical(mars_adata.obs[celltype_col])
    adata_mars.obs.to_csv(result_dir+os.sep+prefix+"_predicted_obs.csv")

    sc.pp.neighbors(adata_mars, n_neighbors=30, use_rep='X')
    sc.tl.umap(adata_mars)
    sc.pl.umap(adata_mars, color=['ground_truth','MARS_labels'],size=15)
    plt.savefig(result_dir+os.sep+prefix+"_prediction_result.png")
    print("=== Finish visualizing..")

def process_data(data_dir, result_dir, args=None):
    '''Load processed datasets from other methods

    I do not want to repeat the configuration generation procedure here for every experiment
    '''
    train_adata, test_adata = None, None
    train_adata_file = result_dir+os.sep+"train_adata.h5ad"
    if os.path.exists(train_adata_file):
        train_adata = anndata.read_h5ad(train_adata_file)

    test_adata_file = result_dir+os.sep+"test_adata.h5ad"
    if os.path.exists(test_adata_file):
        test_adata = anndata.read_h5ad(test_adata_file)

    if train_adata is None or test_adata is None:
        sys.exit("Please run other methods to generate the train/test anndata.")

    ## find common celltypes between train and test
    common_celltypes = set(train_adata.obs["cell.type"]).intersection(set(test_adata.obs["cell.type"]))
    train_cells = train_adata.obs.loc[train_adata.obs["cell.type"].isin(common_celltypes)].index
    test_cells = test_adata.obs.loc[test_adata.obs["cell.type"].isin(common_celltypes)].index

    train_adata = train_adata[train_cells]
    test_adata = test_adata[test_cells]
    return train_adata, test_adata


