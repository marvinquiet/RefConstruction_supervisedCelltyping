import os
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
from textwrap import wrap


def plot_metrics_bar(result_dirs, plot_dir, metrics="Acc", prefix="pancreas",
        groupby="features"):
    '''Plot Acc, ARI, macroF1 for different methods
    x-axis is different feature selection
    y-axis is different matrics'

    @result_dirs: given a list of result_dir
    @metrics: one of evaluation metrics Acc/ARI/macroF1
    @groupby: 
        - features: group by features and compare between differnt feature selections
        - methods: group by methods and compare between methods
    '''
    res_dict = {}
    for result_dir in result_dirs:
        basename = os.path.basename(result_dir)
        basename = basename.replace("result_"+prefix+'_', '')
        res_dict[basename] = {}
        files = [x for x in os.listdir(result_dir) if x.endswith("_metrics.txt")]
        methods = [x.replace("_metrics.txt", "") for x in files] ## get all methods
        for method in methods:
            res_dict[basename][method] = {}
            with open(result_dir+os.sep+method+"_metrics.txt", 'r') as fopen:
                for line in fopen:
                    metric, value = line.split(":")
                    if metric == metrics:
                        res_dict[basename][method] = float(value)
    df = pd.DataFrame.from_dict(res_dict)
    df = df.reindex(sorted(df.columns), axis=1)
    #df.columns = ['\n'.join(wrap(x, 14)) for x in  df.columns] ## cut into shorter

    ## line plot
    #df.T.plot.line(rot=45)

    ## === for all methods
    #level=["scmap", "CHETAH", "RF", "SVM_linear", "SVM_RBF", "MLP", "MLP_GO", "MLP_CP", "DFN", "GEDFN", "ItClust"]
    #df = df.reindex(level)

    fig = plt.figure(figsize=(10, 6), dpi=300)
    ax = plt.gca()
    ## bar plot
    if "features" == groupby:
        if metrics in ["Acc", "ARI", "macroF1"]:
            df.T.plot(ax=ax, kind="bar", rot=45)
        else:
            df.T.plot(ax=ax, kind="bar", rot=45, logy=True)
        plt.xlabel("Feature Selection")
    elif "methods" == groupby:
        if metrics in ["Acc", "ARI", "macroF1"]:
            df.plot(ax=ax, kind="bar", rot=90)
        else:
            df.plot(ax=ax, kind="bar", rot=90, logy=True)
        plt.xlabel("Classifiers")

    ## set y axis range
    global_min = df.to_numpy().min()
    global_max = df.to_numpy().max()
    if "Acc" == metrics or "macroF1" == metrics:
        min_y_axis = max(global_min*0.95, 0)
        max_y_axis = min(global_max*1.05, 1)
        plt.ylim([min_y_axis, max_y_axis])
    elif "ARI" == metrics:
        min_y_axis = max(global_min*0.95, -1)
        max_y_axis = min(global_max*1.05, 1)
        plt.ylim([min_y_axis, max_y_axis])

    plt.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
    plt.ylabel(metrics)
    plt.tight_layout()
    plt.savefig(plot_dir+os.sep+prefix+metrics+'.png')

def extract_sub_prediction(result_dirs, datatype="mousebrain"):
    ''' Extract prediction information from mousebrain sub-cell types

    @datatype: mousebrain/humanPBMC -> belongs to two-stage analysis
               mousebrain_sub/humanPBMC_sub -> belongs to directly prediction using sub cell types
    '''
    if "mousebrain" == datatype:
        major_celltype_col = "mouse_celltypes"
        sub_celltype_col = "cell.type"
        pred_celltype_col = "pred_sub_celltypes"
    elif "mousebrain_sub" == datatype:
        major_celltype_col = "mouse_celltypes"
        sub_celltype_col = "cell.type"
        pred_celltype_col = "pred_celltypes"
    elif "humanPBMC" == datatype:
        major_celltype_col = "cell.type"
        sub_celltype_col = "subtypes"
        pred_celltype_col = "pred_sub_celltypes"
    elif "humanPBMC_sub" == datatype:
        major_celltype_col = "majortypes"
        sub_celltype_col = "cell.type"
        pred_celltype_col = "pred_celltypes"

    import json
    from sklearn import metrics
    for result_dir in result_dirs:
        basename = os.path.basename(result_dir)
        basename = basename.replace("result_"+prefix+'_', '')

        suffix="_predicted_obs.csv"
        files = [x for x in os.listdir(result_dir) if x.endswith(suffix)]
        methods = [x.replace(suffix, "") for x in files] ## get all methods

        for method in methods:
            json_list = []
            df = pd.read_csv(result_dir+os.sep+method+suffix)
            for celltype in set(df[major_celltype_col]):
                json_dict = {}
                json_dict["major_celltypes"] = celltype

                ## for each cell type calculate number of samples accuracy, ARI, macroF1
                sub_df = df[df[major_celltype_col] == celltype]

                ## total major cell numbers
                json_dict["total_number"] = sub_df.shape[0]

                ## sub-cell type numbers
                sub_counts = sub_df[sub_celltype_col].value_counts()
                sub_celltypes_str = ""
                for index in sub_counts.index:
                    sub_celltypes_str += str(index)
                    sub_celltypes_str += ":"
                    sub_celltypes_str += str(sub_counts[index])
                    sub_celltypes_str += "\n"
                json_dict["sub_celltypes"] = sub_celltypes_str

                ## evaluation metrics
                acc = metrics.accuracy_score(sub_df[sub_celltype_col], sub_df[pred_celltype_col])
                ARI = metrics.cluster.adjusted_rand_score(sub_df[sub_celltype_col], sub_df[pred_celltype_col])
                macroF1 = metrics.f1_score(sub_df[sub_celltype_col], sub_df[pred_celltype_col], average="macro")
                json_dict["Accuracy"] = round(acc, 3)
                json_dict["ARI"] = round(ARI, 3)
                json_dict["macroF1"] = round(macroF1, 3)
                json_list.append(json_dict)

            with open(result_dir+os.sep+method+"_summary.json", 'w') as f:
                json.dump(json_list, f)

def plot_feature_number(pipeline_dir, prefix, method="MLP", metrics="Acc"):
    '''Line plot for different feature numbers

    @pipeline_dir: directory of pipeline
    @prefix: experiment prefix
    @method: MLP
    @metrics: Acc/ARI/macroF1
    '''
    res_dict = {}
    for i in range(1, 51):
        n_features = i*100
        result_dir = pipeline_dir+os.sep+'result_'+prefix+'_FEAST_F-test_'+str(n_features)+'_on_train'
        with open(result_dir+os.sep+method+"_metrics.txt", 'r') as fopen:
            for line in fopen:
                metric, value = line.split(":")
                if metric == metrics:
                    res_dict[n_features] = float(value)

    fig = plt.figure(figsize=(10, 6), dpi=300)
    plt.plot(res_dict.keys(), res_dict.values(), 'o', color="red", linestyle='dashed')
 
    plt.xlabel("Number of features")
    plt.ylabel(metrics)
    plt.tight_layout()
    plt.savefig(pipeline_dir+os.sep+prefix+metrics+'_n_features.png')

if __name__ == '__main__':
    pipeline_dir = "/home/wma36/gpu/sc_identifier/pipelines/"
    sub_dirs = next(os.walk(pipeline_dir))[1]

    prefix="PBMC_batch1_ABC_A_to_B"

    result_dirs = [pipeline_dir+os.sep+x for x in sub_dirs if prefix in x]

    groupby = "methods"  ## features/methods
    plot_metrics_bar(result_dirs, pipeline_dir, metrics="Acc", prefix=prefix, groupby=groupby)
    plot_metrics_bar(result_dirs, pipeline_dir, metrics="ARI", prefix=prefix, groupby=groupby)
    plot_metrics_bar(result_dirs, pipeline_dir, metrics="macroF1", prefix=prefix, groupby=groupby)
    plot_metrics_bar(result_dirs, pipeline_dir, metrics="runtime", prefix=prefix, groupby=groupby)

    extract_sub_prediction(result_dirs, datatype="mousebrain_sub")

    method = "MLP"
    #plot_feature_number(pipeline_dir, prefix, method=method, metrics="Acc")
    #plot_feature_number(pipeline_dir, prefix, method=method, metrics="ARI")
    #plot_feature_number(pipeline_dir, prefix, method=method, metrics="macroF1")
    #plot_feature_number(pipeline_dir, prefix, method=method, metrics="runtime")
