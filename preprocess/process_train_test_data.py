import os, sys
import numpy as np
import pandas as pd
import scipy
import scanpy.api as sc
import anndata

## plot
import matplotlib.pyplot as plt

from preprocess import load_pancreatic_data
from preprocess import load_mousebrain_data
from preprocess import load_mouseFC_data
from preprocess import purify_cells

FEAST_SC3_RSCRIPT_PATH = "/homelocal/wma36/celltyping_refConstruct/pipelines/Rcode/FEAST_selection.R"
FEAST_FTEST_RSCRIPT_PATH = "/homelocal/wma36/celltyping_refConstruct/pipelines/Rcode/Ftest_selection.R"

## ---- some functions for processing data
def process_adata(adata, min_genes=10, min_cells=10, celltype_label="cell.type"):
    '''Procedures for filtering single-cell data
       1. Filter low-quality cells and genes;
       2. Filter nonsense genes;
       3. Normalize and log-transform the data;
       4. Change all gene names into UPPER;
       5. Remove cells with no labels; 
    '''
    adata.var_names=[i.upper() for i in list(adata.var_names)]#avod some genes having lower letter

    ## make names unique
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    #1.pre filter cells
    # prefilter_cells(adata,min_genes=min_genes) 
    sc.pp.filter_cells(adata, min_genes=min_genes)

    #2 pre_filter genes
    #prefilter_genes(adata,min_cells=min_cells) # avoiding all gene is zeros
    sc.pp.filter_genes(adata, min_cells=min_genes)

    #3 prefilter_specialgene: MT and ERCC  -> from ItClust package
    Gene1Pattern="ERCC"
    Gene2Pattern="MT-"
    id_tmp1=np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names],dtype=bool)
    id_tmp2=np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names],dtype=bool)
    id_tmp=np.logical_and(id_tmp1,id_tmp2)
    adata._inplace_subset_var(id_tmp)

    ## handel exception when there are not enough cells or genes after filtering
    if adata.shape[0] < 3 or adata.shape[1] < 3:
        return None

    #4 normalization,var.genes,log1p
    sc.pp.normalize_per_cell(adata)  ## total count equal to the median of the counts_per_cell
    sc.pp.log1p(adata)

    ## cells with celltypes
    cells = adata.obs.dropna(subset=[celltype_label]).index.tolist()
    adata = adata[cells]
    return adata

def feature_selection_train_test(train_adata, test_adata, result_dir,
        gene_no=1000, select_on="test", select_method="Seurat",
        min_genes=10, min_cells=10, celltype_label="cell.type"):
    '''Perform feature selection on train_adata and test_adata

    @ result_dir: when method is FEAST, for storing temp counts and FEAST features
    @ gene_no: top x numbers of features
    @ select_on: whether perform on test or train
    @ select_method: Seurat/FEAST/F-test
        FEAST(unsupervised)/F-test(supervised) based on FEAST implementation and
        can be only applied to training datasets because we do not know labels
        for test
    @ min_genes/min_cells: when processing the data, the minimun requirements
    '''

    ## if feature already exists
    feature_file = result_dir+os.sep+"features.txt"
    if os.path.exists(feature_file):
        ## filter cells/genes, etc
        train_adata = process_adata(train_adata, min_genes, min_cells)
        test_adata = process_adata(test_adata, min_genes, min_cells)

        with open(feature_file) as f:
            features = f.read().splitlines()
        print("Number of features:", len(features))

        features = set(train_adata.var_names.tolist()).intersection(set(features))
        features = set(test_adata.var_names.tolist()).intersection(set(features))
        features = list(features)
        features.sort()  ## for reproducibility

        ## order common genes in anndata
        train_adata = train_adata[:, features]
        test_adata = test_adata[:, features]
        return train_adata, test_adata

    if "FEAST" == select_method or "F-test" == select_method:
        ## write file to result dir and run Rscript
        if select_on == "test":
            tmp_adata = test_adata.copy()
        if select_on == "train":
            tmp_adata = train_adata.copy()
        ## to csv
        tmp_data = None
        if scipy.sparse.issparse(tmp_adata.X) or \
                isinstance(tmp_adata.X, pd.DataFrame):
            tmp_data = tmp_adata.X.toarray()
        else:
            tmp_data = tmp_adata.X

        ## write out original read count matrix
        tmp_df = pd.DataFrame(data=tmp_data, index=tmp_adata.obs_names, columns=tmp_adata.var_names).T
        tmp_df_path = result_dir+os.sep+"tmp_counts.csv"
        tmp_df.to_csv(tmp_df_path)

        if "FEAST" == select_method:
            os.system("Rscript --vanilla " + FEAST_SC3_RSCRIPT_PATH + " "+ tmp_df_path + " " + 
                    str(len(set(tmp_adata.obs[celltype_label]))) + " " + str(gene_no))
        elif "F-test" == select_method and "train" == select_on:
            ## write out cell annotations based on train
            cell_annots = tmp_adata.obs[celltype_label].tolist()
            cell_annots_path = result_dir+os.sep+"tmp_cell_annots.txt"
            with open(cell_annots_path, 'w') as f:
                for cell_annot in cell_annots:
                    f.write("%s\n" % cell_annot)
            os.system("Rscript --vanilla " + FEAST_FTEST_RSCRIPT_PATH + " "+ tmp_df_path + " " + 
                    cell_annots_path + " " + str(gene_no))
            os.system("rm {}".format(cell_annots_path))  ## remove the temporaty cell annotations

        os.system("rm {}".format(tmp_df_path))  ## remove the temporaty counts
        del tmp_adata

    ## filter cells/genes, etc
    train_adata = process_adata(train_adata, min_genes, min_cells)
    test_adata = process_adata(test_adata, min_genes, min_cells)

    ## handle None exception
    if train_adata is None or test_adata is None:
        return None, None

    ## select top 1000 HVGs from test
    if select_method == "Seurat":
        if select_on == "test":
            sc.pp.highly_variable_genes(test_adata, n_top_genes=gene_no, subset=True)
        if select_on == "train":
            sc.pp.highly_variable_genes(train_adata, n_top_genes=gene_no, subset=True)

    if "FEAST" == select_method or "F-test" == select_method:
        ## read features selected by FEAST
        feast_file = result_dir+os.sep+select_method+"_features.txt"
        with open(feast_file) as f:
            feast_features = f.read().splitlines()
        feast_features = [x.upper() for x in feast_features] ## upper case

        if select_on == "test":
            feast_genes = set(feast_features).intersection(set(test_adata.var_names.tolist()))
            test_adata = test_adata[:, list(feast_genes)]
        if select_on == "train":
            feast_genes = set(feast_features).intersection(set(train_adata.var_names.tolist()))
            train_adata = train_adata[:, list(feast_genes)]

    features = set(train_adata.var_names.tolist()).intersection(set(test_adata.var_names.tolist()))
    features = list(features)
    features.sort()  ## for reproducibility
    print("Number of features:", len(features))

    ## write features into file
    with open(result_dir+os.sep+"features.txt", 'w') as f:
        for feature in features:
            f.write("%s\n" % feature)

    ## order common genes in anndata
    train_adata = train_adata[:, features]
    test_adata = test_adata[:, features]
    return train_adata, test_adata

def scale_and_visualize(train_adata, test_adata, result_dir, dr_seed=0, scale=True,
        plot=True, plot_elements=['dataset_batch', 'cell.type'], 
        purify_method="", purify_rate=0.1):
    '''Scale data set and plot a dimension reduction on certain elements
    @dr_seed: seed for dimention reduction
    @plot_elements: dataset_batch, cell.type, or ind
    @purify: whether to purify the data or not
    @purify_method: distance or SVM -> whether to use distance or fit probability to purify
    @purify_rate: how many to drop
    '''
    if scale:
        sc.pp.scale(train_adata, zero_center=True, max_value=6)
        sc.pp.scale(test_adata, zero_center=True, max_value=6)

    if purify_method != "":
        if purify_method == "distance":
            train_adata = purify_cells.purify_distances(train_adata, drop_rate=purify_rate)
        elif purify_method == "SVM":
            train_adata = purify_cells.purify_SVM(train_adata, drop_rate=purify_rate)
        plot_adata(train_adata, ['remained'], result_dir, dr_seed=dr_seed, 
                prefix="train_remained_")
        train_adata = train_adata[train_adata.obs["remained"] == True] # use remained data

        ## curate for common cell types in train/test adata
        common_celltypes = set(train_adata.obs["cell.type"]).intersection(set(test_adata.obs["cell.type"]))
        train_cells = train_adata.obs.loc[train_adata.obs["cell.type"].isin(common_celltypes)].index
        test_cells = test_adata.obs.loc[test_adata.obs["cell.type"].isin(common_celltypes)].index
        train_adata = train_adata[train_cells]
        test_adata = test_adata[test_cells]

    ## plot train_adata and test_adata
    adata = train_adata.concatenate(test_adata,join='inner',
            batch_key="dataset_batch",batch_categories=["train","test"]) #inner join

    if plot:
        if "ind" in adata.obs.columns:
            adata.obs["ind"] = adata.obs["ind"].astype("category")
        plot_adata(adata, plot_elements, result_dir, dr_seed=dr_seed)

    ## set adata information to train_adata, test_adata
    train_adata = adata[adata.obs[adata.obs["dataset_batch"] == "train"].index.tolist()]
    test_adata = adata[adata.obs[adata.obs["dataset_batch"] == "test"].index.tolist()]
    return train_adata, test_adata

def load_adata(result_dir):
    '''If data already exists, load from disk
    '''
    if (os.path.exists(result_dir+os.sep+"train_adata.h5ad") and
            os.path.exists(result_dir+os.sep+"test_adata.h5ad")):
        train_adata = anndata.read_h5ad(result_dir+os.sep+"train_adata.h5ad")
        test_adata = anndata.read_h5ad(result_dir+os.sep+"test_adata.h5ad")
        return True, train_adata, test_adata
    else:
        return False, None, None
 
def save_adata(train_adata, test_adata, result_dir, write=True):
    '''Save data to disk
    '''
    if write:
        train_adata.write(result_dir+os.sep+"train_adata.h5ad")
        test_adata.write(result_dir+os.sep+"test_adata.h5ad")

def plot_adata(adata, columns, result_dir, dr_seed=0, prefix=""):
    '''Dimension reduction on adata and plot out features of interest

    @ adata: combined anndata of train and test
    @ columns: columns of interest from adata.obs 
    @ result_dir: where to store the plots
    @ dr_seed: dimension reduction seed
    '''
    # do PCA first
    sc.tl.pca(adata, random_state=dr_seed)
    sc.tl.tsne(adata, use_rep="X_pca",
            learning_rate=300, perplexity=30, n_jobs=4, random_state=dr_seed)
    sc.pl.tsne(adata, color=columns)
    plt.savefig(result_dir+os.sep+prefix+"tSNE_cluster.png")

    ## do UMAP
    sc.pp.neighbors(adata, n_neighbors=20, use_rep="X_pca", random_state=dr_seed) 
    sc.tl.umap(adata, random_state=dr_seed)
    sc.pl.umap(adata, color=columns)
    plt.savefig(result_dir+os.sep+prefix+"umap_cluster.png")

def process_pipeline(train_adata, test_adata, result_dir, 
        gene_no=1000, select_on="test", select_method="Seurat",
        min_genes=10, min_cells=10):
    ''' A process pipeline integrating feature selection, center scaled the data
    '''
    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method, min_genes, min_cells)

    ## if after feature selection, one of them is None, then train and test to None
    if train_adata is None or test_adata is None:
        return None, None

    ## scale and analyze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir,
            plot=False)
    return train_adata, test_adata


