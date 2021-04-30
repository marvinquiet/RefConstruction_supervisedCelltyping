import os,sys
import numpy as np
import pandas as pd
import scipy

import anndata

def purify_distances(adata, drop_rate=0.1):
    '''Purify the cells using distance

    @adata: anndata object
    @drop_rate: drop the top xx distances

    return: original adata, but with one more column in obs dataframe
    '''
    celltype_ind = "cell.type"

    adata_obs = adata.obs
    remained_cells = []
    for cell_type in set(adata.obs[celltype_ind]):
        cells = adata_obs[adata_obs[celltype_ind] == cell_type].index
        sub_adata = adata[cells]
        centroid = np.mean(sub_adata.X, axis=0)  ## row mean, which is the centroid
        distances = np.linalg.norm(sub_adata.X - centroid, axis=1) ## distances for each cell
        top_k_num = round(len(distances)*drop_rate)
        if 0 == top_k_num: ## if no need to remove, continue
            continue

        top_k_indices = distances.argsort()[-top_k_num:][::-1] ## get top indexes
        remained_index = sub_adata.obs.index.delete(top_k_indices.tolist())
        remained_cells.extend(remained_index.tolist())

    adata.obs["remained"] = [True if cell in remained_cells else False for cell in adata.obs.index]
    adata.obs["remained"] = adata.obs["remained"].astype("category")
    return adata

def purify_SVM(adata, drop_rate=0.1, kernel="linear", seed=0):
    '''Use SVM to decide probability to drop
    
    Note: to speed up, use linear kernel
    '''
    celltype_ind = "cell.type"

    from sklearn.svm import SVC
    ## apply 5-fold validation to get probabilities
    model = SVC(decision_function_shape="ovr", kernel=kernel, random_state=seed,
            probability=True)

    from sklearn.preprocessing import OneHotEncoder

    enc = OneHotEncoder(handle_unknown='ignore')
    x_train = np.array(adata.X)
    y_train = enc.fit_transform(adata.obs[[celltype_ind]]).toarray()

    model.fit(x_train, y_train.argmax(1))
    y_pred = model.predict_proba(x_train)
    labels = y_train.argmax(1)
    celltypes = adata.obs[celltype_ind].tolist()
    mapping = dict(zip(celltypes, labels))

    remained_cells = []
    for cell_type in set(adata.obs[celltype_ind]):
        idx = np.where(adata.obs[celltype_ind] == cell_type)[0] ## cells matched with the cell type
        cell_pred = y_pred[idx, mapping[cell_type]]

        top_k_num = round(len(cell_pred)*drop_rate)
        top_k_indices = cell_pred.argsort()[: top_k_num] ## get the most not confident cells
        remained_index = np.delete(idx, top_k_indices)
        remained_cells.extend(remained_index)
    adata.obs["remained"] = [True if i in remained_cells else False for i in range(len(adata.obs.index))]
    adata.obs["remained"] = adata.obs["remained"].astype("category")
    return adata
