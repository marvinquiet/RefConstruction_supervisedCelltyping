'''
A celltyping pipeline integrating DFN, GEDFN, MLP, SVM and Random Forest

For scmap and CHETA, there are other Rscripts and then use bash to run
'''

import os, sys, math, time
import anndata
import numpy as np
import pandas as pd
import scipy
import tensorflow as tf

from sklearn.preprocessing import OneHotEncoder
from sklearn import metrics

from random import seed
from numpy.random import seed
from tensorflow import set_random_seed

# GPU settings and reproducible
#os.environ["CUDA_VISIBLE_DEVICES"]="0"

#Set seeds
RANDOM_SEED=0
seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
set_random_seed(RANDOM_SEED)

## every time reset all graphs and run something new
tf.compat.v1.reset_default_graph()

def run_GEDFN(x_train, y_train, x_test, partition, 
        dims=[64], batch_size=128, GE=True):
    '''Pipeline for running GEDFN

    @ train_adata/test_adata: train/test anndata
    @ partition: gene network
    @ dims: for GEDFN and DFN, the dimension of the second to the last (because the
    1st hidden layer has the same dimension as input)
    @ celltype_cols: the column indicator of cell type
    '''
    from test_GEDFN.test_GEDFN import GEDFN
    ## initiate GEDFN
    gedfn = GEDFN(dims=dims, partition=partition, GE=GE)

    ## open tf session for running
    sess = tf.Session()
    gedfn.train(x_train, y_train, tf_session=sess, batch_size=batch_size)

    ## print weights between input and 1st hidden
    #print(sess.run(gedfn.weights['h1']))

    ## prediction and evaluation
    y_pred = gedfn.predict(x_test, tf_session=sess)
    # print("=== Predicted Y:", y_pred)

    ## extract important features
    print("=== Extracting important features...")
    var_importance_score = gedfn.feature_importance(sess)
    # print("=== Feature importance: ", var_important)

    ## close tf session
    sess.close()
    return y_pred, var_importance_score

def run_MLP(x_train, y_train, x_test, dims=[128, 64, 32, 16], 
        batch_size=128, seed=0):
    '''Pipeline for MLP
    '''

    from test_GEDFN.test_MLP import MLP

    mlp = MLP(dims)
    mlp.fit(x_train, y_train, batch_size=batch_size, random_state=seed)
    y_pred = mlp.predict(x_test)
    return y_pred

def run_GEMLP(x_train, y_train, x_test, partition,
        dims=[128, 64, 32, 16], batch_size=128):
    from test_GEDFN.test_GEMLP import GEMLP

    ## initialize GEMLP
    gemlp = GEMLP(dims=dims, partition=partition)
    gemlp.fit(x_train, y_train, batch_size=batch_size)
    y_pred = gemlp.predict(x_test)
    return y_pred


def run_ItClust(x_train, y_train, x_test, 
        dims=[64], batch_size=128, seed=0):
    '''Pipeline for ItClust
    '''
    from ItClust.DEC import DEC
    from keras.optimizers import SGD
    ## parameters
    ALPHA = 1.0
    INIT = "glorot_uniform"
    PRETRAIN_EPOCHS = 300
    TOL = 0.001

    dims.insert(0, x_train.shape[1])

    ## first round DEC
    dec=DEC(dims=dims, y=pd.Series(y_train.argmax(1)), x=x_train,
            alpha=ALPHA,init=INIT,pretrain_epochs=PRETRAIN_EPOCHS,
            actinlayer1="tanh",softmax=False)
    dec.compile(optimizer=SGD(lr=0.01,momentum=0.9))
    Embeded_z,q_pred=dec.fit_supervise(x = x_train,y = y_train.argmax(1), epochs=2e3,batch_size=batch_size)
    ## transfer weights
    weights=[i0.get_weights() for i0 in dec.model.layers]
    features=dec.encoder.predict(x_test)
    q=dec.model.predict(x_test,verbose=0)
    ## check training result
    print("=== Training model finished! Start to fit target network!")
    val_y_pre=dec.model.predict(x_train,verbose=0)
    val_y_pre=[np.argmax(i) for i in val_y_pre]
    ARI=metrics.adjusted_rand_score(y_train.argmax(1), val_y_pre)
    acc = metrics.accuracy_score(y_train.argmax(1), val_y_pre)
    macroF1 = metrics.accuracy_score(y_train.argmax(1), val_y_pre)
    print("*****=====", "Training accuracy: ", acc, " Training ARI: ", ARI, 
            "Trainig MacroF1: ", macroF1, "=====*****")
    ## fit target
    dec2=DEC(dims=dims,x=x_test,
            alpha=ALPHA,init=INIT,pretrain_epochs=PRETRAIN_EPOCHS,
            actinlayer1="tanh",softmax=False,transfer_feature=features,
            model_weights=weights,y_trans=q.argmax(axis=1))
    dec2.compile(optimizer=SGD(0.01,0.9))
    trajectory_z, trajectory_l, Embeded_z,q_pred=dec2.fit_trajectory(x=x_test,
            tol=TOL,epochs_fit=5,batch_size=batch_size)# Fine tunning

    return q_pred


def run_SVM(x_train, y_train, x_test, kernel="rbf", seed=0):
    '''Fit SVM model with a RBF kernel
    SVM decision_function_shape: can be one-vs-one or one-vs-rest, in our case,
    according to our input, we should use one-vs-rest
    '''
    if "rbf" == kernel:
        from sklearn.svm import SVC
        model = SVC(decision_function_shape="ovr", kernel=kernel, random_state=seed)
        #model = SVC(decision_function_shape="ovo", kernel=kernel, random_state=seed)
    elif "linear" == kernel:
        from sklearn.svm import LinearSVC
        model = LinearSVC(multi_class='ovr', max_iter=1e4, random_state=seed)

    ## fit model
    model.fit(x_train, y_train)
    y_pred = model.predict(x_test)
    return y_pred

def run_RF(x_train, y_train, x_test, tree=50, seed=0):
    '''Fit Random Forest
    '''
    from sklearn.ensemble import RandomForestClassifier
    model = RandomForestClassifier(n_estimators=50, random_state=seed)

    ## fit model
    model.fit(x_train, y_train)
    y_pred = model.predict(x_test)
    return y_pred

def analyze_prediction(y_pred, y_test, enc_obj, test_adata, n_clusters, 
        run_time=None, celltype_cols="cell.type", result_dir="./result", prefix="model"):
    ''' Analyze the celltyping prediction result
    @ y_pred: the result from predicting
    @ enc_obj: onehot encoder object to convert back to cell types
    @ n_clusters: number of clusters in train datasets
    @ thres: if thres, then we generate "unassigned" results
    '''
    ## evaluation metrics
    acc = metrics.accuracy_score(y_test, y_pred)
    ARI = metrics.cluster.adjusted_rand_score(y_test, y_pred)
    macroF1 = metrics.f1_score(y_test, y_pred, average="macro")
    print("*****=====", "Testing accuracy: ", acc, " Testing ARI: ", ARI, 
            "Testing MacroF1: ", macroF1, "=====*****")
    
    ## write metrics result to file
    with open(result_dir+os.sep+prefix+"_metrics.txt", 'w') as f:
        f.write("Acc:%s\n" % str(acc))
        f.write("ARI:%s\n" % str(ARI))
        f.write("macroF1:%s\n" % str(macroF1))
        if run_time is not None:
            f.write("runtime:%s\n" % str(run_time))

    ## analyze predicted cell types 
    pred_labels = y_pred
    pred_onehot = np.zeros((pred_labels.size, n_clusters))
    pred_onehot[np.arange(pred_labels.size), pred_labels] = 1
    pred_celltypes = enc_obj.inverse_transform(pred_onehot)
    print("=== Predicted celltypes: ", set(pred_celltypes.flatten()))

    test_adata.obs['pred_celltypes'] = pred_celltypes
    test_adata.obs.to_csv(result_dir+os.sep+prefix+"_predicted_obs.csv")
 
    ## visualization
    import matplotlib.pyplot as plt
    import scanpy.api as sc
    sc.pl.tsne(test_adata, color=[celltype_cols, "pred_celltypes"], size=15)
    plt.savefig(result_dir+os.sep+prefix+"_prediction_result.png")
    print("=== Finish visualizing..")

def run_pipeline(args, train_adata, test_adata, data_dir, result_dir):
    '''Run methods
    '''
    ## Hyperparameters for network
    if train_adata.shape[0] >= 5000:
        ## === parameters for mousebrain (high cell number)
        dims = [128, 32]
        MLP_dims = [128, 64, 32, 16, 8]
    else:
        ## === parameters for PBMC datasets (low cell number)
        dims = [16]
        MLP_dims = [64, 16]

    batch_size = 128
    celltype_cols = "cell.type"
    n_clusters = len(set(train_adata.obs[celltype_cols]))
    ## OneHotEncoding the celltypes
    enc = OneHotEncoder(handle_unknown='ignore')
    if scipy.sparse.issparse(train_adata.X):
        x_train = train_adata.X.toarray()
    else:
        x_train = train_adata.X
    y_train = enc.fit_transform(train_adata.obs[[celltype_cols]]).toarray()

    if scipy.sparse.issparse(test_adata.X):
        x_test = test_adata.X.toarray()
    else:
        x_test = test_adata.X
    y_test = enc.transform(test_adata.obs[[celltype_cols]]).toarray()

    if "MLP" == args.method:
        ### --- run MLP
        print("\n\n=== MLP\n")
        start = time.time()
        y_pred = run_MLP(x_train, y_train, x_test, 
                dims=MLP_dims, batch_size=batch_size, seed=RANDOM_SEED) ## run MLP
        end = time.time()
        print("\n\n=== Run time:", end-start)
        analyze_prediction(y_pred.argmax(1), y_test.argmax(1), enc, test_adata, 
                n_clusters, end-start, celltype_cols=celltype_cols, 
                result_dir=result_dir, prefix=args.method)

    if "MLP_GO" == args.method or "MLP_CP" == args.method:
        ### --- run GEMLP
        network = args.method.split('_')[1]

        if "GO" == network:
            embedded_network = pd.read_csv(data_dir+os.sep+"GeneNetworks/c5.go.v7.2.symbolsadj.txt", sep=" ", index_col=0)
        elif "CP" == network:
            embedded_network = pd.read_csv(data_dir+os.sep+"GeneNetworks/c2.cp.v7.2.symbolsadj.txt", sep=" ", index_col=0)

        embedded_network.index = embedded_network.index.str.upper()
        common_genes = set(embedded_network.index).intersection(set(train_adata.var_names))
        common_genes = list(common_genes)
        common_genes.sort()
    
        partition = embedded_network.loc[common_genes, :]
        partition = partition.to_numpy()  ## turn dataframe to numpy
    
        ## input intersect with partition
        x_train = np.array(train_adata[:, common_genes].X)
        x_test = np.array(test_adata[:, common_genes].X)
        print("After intersecting with partition: ", x_train.shape, x_test.shape)
        y_train = enc.fit_transform(train_adata.obs[[celltype_cols]]).toarray()
        y_test = enc.transform(test_adata.obs[[celltype_cols]]).toarray()

        print("\n\n=== %s\n" % args.method)
        start = time.time()
        y_pred = run_GEMLP(x_train, y_train, x_test, partition,
                dims=MLP_dims, batch_size=batch_size) ## run GEMLP
        end = time.time()
        print("\n\n=== Run time:", end-start)
        analyze_prediction(y_pred.argmax(1), y_test.argmax(1), enc, test_adata, 
                n_clusters, end-start, celltype_cols=celltype_cols, 
                result_dir=result_dir, prefix=args.method)
 

    if "GEDFN" == args.method or "DFN" == args.method:
        ## load partition dataframe
        partition_df = pd.read_csv(data_dir+os.sep+"GeneNetworks/HomoSapients_htb_hq_adj.txt", 
                sep=" ", index_col=0)
        partition_df.index = partition_df.index.str.upper()
        partition_df.columns = partition_df.columns.str.upper()

        ## read features
        with open(result_dir+os.sep+"features.txt", 'r') as f:
            features = f.read().splitlines()

        common_genes = set(features).intersection(set(partition_df.index))
        common_genes = list(common_genes)
        common_genes.sort() ## for reproducibility
        partition = partition_df.loc[common_genes, common_genes]
        partition = partition.to_numpy()  ## turn dataframe to numpy

        ## input intersect with partition
        x_train = np.array(train_adata[:, common_genes].X)
        x_test = np.array(test_adata[:, common_genes].X)
        print("After intersecting with partition: ", x_train.shape, x_test.shape)
        y_train = enc.fit_transform(train_adata.obs[[celltype_cols]]).toarray()
        y_test = enc.transform(test_adata.obs[[celltype_cols]]).toarray()

        ### --- run GEDFN
        print("\n\n=== "+args.method+"\n")
        start = time.time()
        if "GE" in args.method:
            y_pred, var_importance = run_GEDFN(x_train, y_train, x_test, partition, 
                                     dims=dims, batch_size=batch_size, GE=True) ## run GEDFN
        else:
            y_pred, var_importance = run_GEDFN(x_train, y_train, x_test, partition, 
                                     dims=dims, batch_size=batch_size, GE=False) ## run DFN
        end = time.time()
        print("\n\n=== Run time:", end-start)
        analyze_prediction(y_pred.argmax(1), y_test.argmax(1), enc, test_adata, 
                n_clusters, end-start, 
                celltype_cols=celltype_cols, result_dir=result_dir, prefix=args.method)

    if "ItClust" == args.method:
        ## --- run ItClust
        print("\n\n=== ItClust\n")
        start = time.time()
        y_pred = run_ItClust(x_train, y_train, x_test, 
                            dims=dims, batch_size=batch_size) ## runItClust
        end = time.time()
        print("\n\n=== Run time:", end-start)
        analyze_prediction(y_pred.argmax(1), y_test.argmax(1), enc, test_adata, 
                n_clusters, end-start, 
                celltype_cols=celltype_cols, result_dir=result_dir, prefix=args.method)
 
    if "SVM_RBF" == args.method:
        ## --- run SVM
        print("\n\n=== SVM RBF kernel\n")
        start = time.time()
        y_pred = run_SVM(x_train, y_train.argmax(1), x_test)
        end = time.time()
        print("\n\n=== Run time:", end-start)
        analyze_prediction(y_pred, y_test.argmax(1), enc, test_adata, 
                n_clusters, end-start, celltype_cols=celltype_cols, 
                result_dir=result_dir, prefix=args.method)

    if "SVM_linear" == args.method:
        ## --- run SVM
        print("\n\n=== SVM linear kernel\n")
        start = time.time()
        y_pred = run_SVM(x_train, y_train.argmax(1), x_test, kernel="linear")
        end = time.time()
        print("\n\n=== Run time:", end-start)
        analyze_prediction(y_pred, y_test.argmax(1), enc, test_adata, 
                n_clusters, end-start, celltype_cols=celltype_cols, 
                result_dir=result_dir, prefix=args.method)

    if "RF" == args.method:
        ## --- run RF
        print("\n\n=== Random Forest\n")
        start = time.time()
        y_pred = run_RF(x_train, y_train, x_test)
        end = time.time()
        print("\n\n=== Run time:", end-start)
        analyze_prediction(y_pred.argmax(1), y_test.argmax(1), enc, test_adata, 
                n_clusters, end-start, 
                celltype_cols=celltype_cols, result_dir=result_dir, prefix=args.method)

