'''
A celltyping pipeline integrating DFN, GEDFN, MLP, SVM and Random Forest

For scmap and CHETA, there are other Rscripts and then use bash to run
'''

import os, sys, math, time
import anndata
import numpy as np
import pandas as pd
import tensorflow as tf
#tf.compat.v1.enable_eager_execution(
#    config=None,
#    device_policy=None,
#    execution_mode=None
#)

from sklearn.utils import shuffle
from sklearn.preprocessing import OneHotEncoder
from sklearn import metrics

from random import seed
from numpy.random import seed
from tensorflow import set_random_seed

# GPU settings and reproducible
#os.environ["CUDA_VISIBLE_DEVICES"]="0"
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


def analyze_prediction(y_pred, y_test, enc_obj, test_adata, n_clusters, run_time=None,
        celltype_cols="cell.type", result_dir="./result", prefix="model"):
    ''' Analyze the celltyping prediction result
    @ y_pred: the result from predicting
    @ enc_obj: onehot encoder object to convert back to cell types
    @ n_clusters: number of clusters in train datasets
    '''
    ## extract number of cell types from test anndata
    #n_clusters = len(set(test_adata.obs[celltype_cols]))

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
    pred_celltypes = enc.inverse_transform(pred_onehot)
    print("=== Predicted celltypes: ", set(pred_celltypes.flatten()))

    test_adata.obs['pred_celltypes'] = pred_celltypes
    test_adata.obs.to_csv(result_dir+os.sep+prefix+"_predicted_obs.csv")
 
    ## visualization
    import matplotlib.pyplot as plt
    import scanpy.api as sc
    sc.pl.tsne(test_adata, color=[celltype_cols, "pred_celltypes"], size=15)
    plt.savefig(result_dir+os.sep+prefix+"_prediction_result.png")
    print("=== Finish visualizing..")


if __name__ == "__main__":
    data_dir = "/home/wma36/gpu/data"

    ## parse arguments
    import argparse
    parser = argparse.ArgumentParser(description="Celltyping pipeline.")
    parser.add_argument('data_source', help="Load which dataset",
        choices=['PBMC_batch1_ind', 'PBMC_batch1_ind_major',
            'PBMC_batch1_ind_purify_dist', 'PBMC_batch1_ind_purify_SVM',
            'PBMC_batch1_ABC', "PBMC_batch1_ABC_purify_dist", "PBMC_batch1_ABC_purify_SVM",
            'PBMC_batch1_batchtoind', 'PBMC_batch1_multiinds', 'PBMC_batch1_multiinds_sample',
            'PBMC_batch2', "PBMC_batch2_purify_dist", "PBMC_batch2_purify_SVM",
            'PBMC_batch1_batch2_ind', 'PBMC_batch1_batch2_ind_major',
            'PBMC_protocols_pbmc1', 'PBMC_protocols_batch_smart',
            'pancreas', 'pancreas_seg_cond', 'pancreas_custom', 
            'pancreas_seg_mix', 'pancreas_multi_to_multi', 
            "mousebrain_FC", "mousebrain_FC_curate", "mousebrain_FC_sub", "mousebrain_FC_sub_purify_dist",
            "mousebrain_FC_sub_purify_SVM", "mousebrain_HC", "mousebrain_HC_sub",
            "mousebrain_FC_multiinds", "mousebrain_FC_multiinds_sub",
            "mousebrain_FC_multiinds_sample", "mousebrain_FC_multiinds_sub_sample",
            "mousebrain_region", "mousebrain_region_sub", 
            "mousebrain_FC_stage", "mousebrain_FC_stage_sub", "mousebrain_FC_datasets",
            "mousebrain_FC_datasets_multiinds", "mousebrain_FC_datasets_multiinds_sample",
            "mousecortex_protocols_mouseFC_ind", 
            "mousebrain_combined_mouseFC"])

    parser.add_argument('-m', '--method', help="Run which method",
        choices=['MLP', 'GEDFN', 'ItClust', 'SVM_RBF', 'SVM_linear', 'RF'], ## remove DFN
        required=True)
    parser.add_argument('--select_on', help="Feature selection on train or test, or None of them",
        choices=['train', 'test'])
    parser.add_argument('--select_method', help="Feature selection method, Seurat/FEAST or None",
            choices=['Seurat', 'FEAST', 'F-test'])
    parser.add_argument('--n_features', help="Number of features selected",
            default=1000, type=int)
    parser.add_argument('--train', help="Specify which as train", required=True)
    parser.add_argument('--test', help="Specify which as test", required=True)
    parser.add_argument('--sample_seed', help="Downsample seed in combined individual effect", 
            default=0, type=int)

    args = parser.parse_args()

    pipeline_dir = "/home/wma36/gpu/sc_identifier/pipelines"
    result_prefix = pipeline_dir+os.sep+"result_"+args.data_source+'_'+\
        args.train+'_to_'+args.test
    if 0 != args.sample_seed:  ## if not 0 for downsampling seed
        result_prefix += '_'
        result_prefix += str(args.sample_seed)
    os.makedirs(result_prefix, exist_ok=True)

    ## create file directory 
    if args.select_on is None and args.select_method is None:
        result_dir = result_prefix+os.sep+"no_feature"
    else:
        result_dir = result_prefix+os.sep+args.select_method+'_'+\
                str(args.n_features)+'_on_'+args.select_on
    os.makedirs(result_dir, exist_ok=True)

    ## === PBMC datasets
    from preprocess import process_PBMC_train_test
    if args.data_source == "PBMC_batch1_ind":
        train_adata, test_adata = \
            process_PBMC_train_test.process_batch1_ind(data_dir, result_dir, 
                ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "PBMC_batch1_ind_major":
        train_adata, test_adata = \
            process_PBMC_train_test.process_batch1_ind(data_dir, result_dir, 
                ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                celltype_gran=0)
    if args.data_source == "PBMC_batch1_ind_purify_dist":
        train_adata, test_adata = \
            process_PBMC_train_test.process_batch1_ind(data_dir, result_dir, 
                ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                purify=True, purify_method="distance")
    if args.data_source == "PBMC_batch1_ind_purify_SVM":
        train_adata, test_adata = \
            process_PBMC_train_test.process_batch1_ind(data_dir, result_dir, 
                ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                purify=True, purify_method="SVM")
    if args.data_source == "PBMC_batch1_ABC":
        train_adata, test_adata = \
            process_PBMC_train_test.process_batch1_ABC(data_dir, result_dir, 
                batch1=args.train, batch2=args.test,
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "PBMC_batch1_ABC_purify_dist":
        train_adata, test_adata = \
            process_PBMC_train_test.process_batch1_ABC(data_dir, result_dir, 
                batch1=args.train, batch2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                purify=True, purify_method="distance")
    if args.data_source == "PBMC_batch1_ABC_purify_SVM":
        train_adata, test_adata = \
            process_PBMC_train_test.process_batch1_ABC(data_dir, result_dir, 
                batch1=args.train, batch2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                purify=True, purify_method="SVM")
    if args.data_source == "PBMC_batch1_batchtoind":
        train_adata, test_adata = \
                process_PBMC_train_test.process_batch1_batchtoind(data_dir, result_dir,
                        batch=args.train, ind=args.test,
                        select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "PBMC_batch1_multiinds":
        train_adata, test_adata = \
                process_PBMC_train_test.process_batch1_multiinds(data_dir, result_dir,
                        sample_ind=args.train, pred_ind=args.test,
                        select_on=args.select_on, select_method=args.select_method,
                        sample=False)
    if args.data_source == "PBMC_batch1_multiinds_sample":
        train_adata, test_adata = \
                process_PBMC_train_test.process_batch1_multiinds(data_dir, result_dir,
                        sample_ind=args.train, pred_ind=args.test,
                        select_on=args.select_on, select_method=args.select_method,
                        sample=True, sample_seed=args.sample_seed)
    if args.data_source == "PBMC_batch2":
        train_adata, test_adata = \
            process_PBMC_train_test.process_batch2_ctrl_stim(data_dir, result_dir, 
                cond1=args.train, cond2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                gene_no=args.n_features)
    if args.data_source == "PBMC_batch2_purify_dist":
        train_adata, test_adata = \
            process_PBMC_train_test.process_batch2_ctrl_stim(data_dir, result_dir, 
                cond1=args.train, cond2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                purify=True, purify_method="distance")
    if args.data_source == "PBMC_batch2_purify_SVM":
        train_adata, test_adata = \
            process_PBMC_train_test.process_batch2_ctrl_stim(data_dir, result_dir, 
                cond1=args.train, cond2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                purify=True, purify_method="SVM")
    if args.data_source == "PBMC_batch1_batch2_ind":
        train_adata, test_adata = \
                process_PBMC_train_test.process_batch1_batch2_ind(data_dir, result_dir,
                        input1=args.train, input2=args.test,
                        select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "PBMC_batch1_batch2_ind_major":
        train_adata, test_adata = \
                process_PBMC_train_test.process_batch1_batch2_ind(data_dir, result_dir,
                        input1=args.train, input2=args.test,
                        select_on=args.select_on, select_method=args.select_method,
                        celltype_gran=0)


    ## === PBMC 7 protocols datasets
    if args.data_source == "PBMC_protocols_pbmc1":
        train_adata, test_adata = \
            process_PBMC_train_test.process_PBMC_protocols_type(data_dir, result_dir, 
                protocols1=args.train, protocols2=args.test, exp="pbmc1",
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "PBMC_protocols_batch_smart":
        train_adata, test_adata = \
            process_PBMC_train_test.process_PBMC_protocols_batch(data_dir, result_dir, 
                batch1=args.train, batch2=args.test, protocol="Smart-seq2",
                select_on=args.select_on, select_method=args.select_method)
 
    from preprocess import process_Pancreas_train_test
    ## === pancrease datasets
    if args.data_source == "pancreas":
        train_adata, test_adata = \
            process_Pancreas_train_test.process_pancreas(data_dir, result_dir, 
                dataset1=args.train, dataset2=args.test, 
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "pancreas_seg_cond":
        train_adata, test_adata = \
            process_Pancreas_train_test.process_pancreas_seg_cond(data_dir, result_dir, 
                cond1=args.train, cond2=args.test, 
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "pancreas_custom":
        train_adata, test_adata = \
            process_Pancreas_train_test.process_pancreas_custom(data_dir, result_dir, 
                seg_cond=args.train, dataset=args.test, 
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "pancreas_seg_mix":
        train_adata, test_adata = \
            process_Pancreas_train_test.process_pancreas_seg_mix(data_dir, result_dir, 
                main_cond=args.train, pred_cond=args.test, 
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "pancreas_multi_to_multi":
        train_adata, test_adata = \
            process_Pancreas_train_test.process_pancreas_multi_to_multi(data_dir, result_dir, 
                cond1=args.train, cond2=args.test,
                select_on=args.select_on, select_method=args.select_method)

    from preprocess import process_mousebrain_train_test
    ## === mouse adult brain
    if args.data_source == "mousebrain_FC":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_data(data_dir, result_dir, 
                region="FC", ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                gene_no=args.n_features)
    if args.data_source == "mousebrain_FC_curate":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_data(data_dir, result_dir, 
                region="FC", ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                gene_no=args.n_features, curate=True)
    if args.data_source == "mousebrain_FC_sub":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_data(data_dir, result_dir, 
                region="FC", ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                celltype_gran=1)
    if args.data_source == "mousebrain_FC_sub_purify_dist":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_data(data_dir, result_dir, 
                region="FC", ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                celltype_gran=1, purify=True, purify_method="distance")
    if args.data_source == "mousebrain_FC_sub_purify_SVM":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_data(data_dir, result_dir, 
                region="FC", ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                celltype_gran=1, purify=True, purify_method="SVM")
    if args.data_source == "mousebrain_FC_multiinds":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_multiinds_data(
                data_dir, result_dir, region="FC", 
                sample_ind=args.train, pred_ind=args.test,
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "mousebrain_FC_multiinds_sub":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_multiinds_data(
                data_dir, result_dir, region="FC", 
                sample_ind=args.train, pred_ind=args.test,
                select_on=args.select_on, select_method=args.select_method,
                celltype_gran=1)
    if args.data_source == "mousebrain_FC_multiinds_sample":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_multiinds_data(
                data_dir, result_dir, region="FC", sample_seed=args.sample_seed,
                sample_ind=args.train, pred_ind=args.test,
                select_on=args.select_on, select_method=args.select_method,
                sample=True)
    if args.data_source == "mousebrain_FC_multiinds_sub_sample":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_multiinds_data(
                data_dir, result_dir, region="FC", sample_seed=args.sample_seed,
                sample_ind=args.train, pred_ind=args.test,
                select_on=args.select_on, select_method=args.select_method,
                celltype_gran=1, sample=True)
    if args.data_source == "mousebrain_HC":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_data(data_dir, result_dir, 
                region="HC", ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "mousebrain_HC_sub":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_data(data_dir, result_dir, 
                region="HC", ind1=args.train, ind2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                celltype_gran=1)
    if args.data_source == "mousebrain_region":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_regions(data_dir, result_dir, 
                region1=args.train, region2=args.test,
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "mousebrain_region_sub":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_regions(data_dir, result_dir, 
                region1=args.train, region2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                celltype_gran=1)

    ## === mousebrain FC
    if args.data_source == "mousebrain_FC_stage":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_FC_stage(data_dir, result_dir, 
                stage1=args.train, stage2=args.test,
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "mousebrain_FC_stage_sub":
        ## sub-celltypes
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_FC_stage(data_dir, result_dir, 
                stage1=args.train, stage2=args.test,
                select_on=args.select_on, select_method=args.select_method,
                celltype_gran=1)
    if args.data_source == "mousebrain_FC_datasets":
        ## predict using dataset1 to dataset2
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_FC_datasets(data_dir, result_dir, 
                dataset1_input=args.train, dataset2_input=args.test,
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "mousebrain_FC_datasets_multiinds":
        ## predict using dataset1 to dataset2
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_FC_datasets(data_dir, result_dir, 
                dataset1_input=args.train, dataset2_input=args.test,
                select_on=args.select_on, select_method=args.select_method)
    if args.data_source == "mousebrain_FC_datasets_multiinds_sample":
        ## predict using dataset1 to dataset2
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousebrain_FC_datasets(data_dir, result_dir, 
                dataset1_input=args.train, dataset2_input=args.test,
                select_on=args.select_on, select_method=args.select_method,
                sample=True, sample_seed=args.sample_seed)

    ## === Mouse cortex data
    if args.data_source == "mousecortex_protocols_mouseFC_ind":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_mousecortex_protocols_mouseFC_data(data_dir, result_dir,
                input1=args.train, input2=args.test, 
                select_on=args.select_on, select_method=args.select_method,
                curate=True)

    ## === Combine mouse brian data to predict FC
    if args.data_source == "mousebrain_combined_mouseFC":
        train_adata, test_adata = \
            process_mousebrain_train_test.process_FCdatasets_mouseFC_data(data_dir, result_dir,
                datasets_ref=args.train, mousebrain_tgt=args.test,
                select_on=args.select_on, select_method=args.select_method,
                curate=True)

    print("Train anndata: \n", train_adata)
    print("Test anndata: \n", test_adata)

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
    x_train = np.array(train_adata.X)
    y_train = enc.fit_transform(train_adata.obs[[celltype_cols]]).toarray()

    x_test = np.array(test_adata.X)
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

    if "GEDFN" == args.method:
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
        y_pred, var_importance = run_GEDFN(x_train, y_train, x_test, partition, 
                                 dims=dims, batch_size=batch_size, GE=True) ## run GEDFN
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
