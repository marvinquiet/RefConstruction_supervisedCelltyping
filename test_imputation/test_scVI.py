import os, sys
import scvi
import anndata
import scanpy.api as sc
import numpy as np

if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Please make sure the result directory is given")

    result_dir = sys.argv[1]
    train_adata = anndata.read_h5ad(result_dir+os.sep+'train_adata.h5ad')
    test_adata = anndata.read_h5ad(result_dir+os.sep+'test_adata.h5ad')

    ## use scVI to impute data
    scvi.data.setup_anndata(train_adata, layer="counts")
    train_model = scvi.model.SCVI(train_adata)
    train_model.train()
    train_denoised = train_model.get_normalized_expression(train_adata, library_size=1e4)
    train_adata.X = np.array(train_denoised)
    sc.pp.scale(train_adata, zero_center=True, max_value=6)

    scvi.data.setup_anndata(test_adata, layer="counts")
    test_model = scvi.model.SCVI(test_adata)
    test_model.train()
    test_denoised = test_model.get_normalized_expression(test_adata, library_size=1e4)
    test_adata.X = np.array(test_denoised)
    sc.pp.scale(test_adata, zero_center=True, max_value=6)

    ## save to the file
    train_adata.write(result_dir+os.sep+'scVI_train_adata.h5ad')
    test_adata.write(result_dir+os.sep+'scVI_test_adata.h5ad')

