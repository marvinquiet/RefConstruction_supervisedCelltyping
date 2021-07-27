import os
import numpy as np
import pandas as pd
from scipy.stats.mstats import mquantiles
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams.update({'font.size': 18})
from textwrap import wrap

def plot_saturation_line(result_dirs, plot_dir, sep_index=0, prefix="saturation",
        columns=None):
    '''Plot Acc, ARI, macroF1 for testing performance saturation

    @result_dirs: given a list of result_dir
    @columns: the x-axis of the line plot
    '''
    metrics = ["Acc", "ARI", "macroF1"]

    res_dict = {}
    for result_dir in result_dirs:
        basename = os.path.basename(result_dir)
        ind_no = len(basename.split('_')) - sep_index ## count for individual number
        res_dict[ind_no] = {}
        ## find metrics file
        for root, dirs, files in os.walk(result_dir):
            for file in files:
                if "MLP_metrics.txt" in file:
                    filepath = os.path.join(root, file)
                    with open(filepath, 'r') as fopen:
                        for line in fopen:
                            metric, value = line.strip().split(":")
                            if metric in metrics:
                                res_dict[ind_no][metric] = float(value)

    df = pd.DataFrame.from_dict(res_dict)
    df = df.reindex(sorted(df.columns), axis=1)

    if columns:
        df.columns = np.cumsum(columns)

    ## line plot on performance
    fig = plt.figure(figsize=(10, 6), dpi=300)
    ax = plt.gca()
    df.T.plot.line(ax=ax)
    plt.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
    if columns:
        plt.xlabel("No. of cells")
    else:
        plt.xlabel("No. of individuals")
    plt.ylabel("Performance")
    plt.tight_layout()
    plt.savefig(plot_dir+os.sep+prefix+'.png')

def plot_downsample_saturation_line(result_dir, prefix="saturation_downsample",
        add_type="cells"):
    '''Plot Acc, ARI, macroF1 for performance saturation on downsampling

    @result_dir: result directory for downsample saturation
    '''
    metrics = ["Acc", "ARI", "macroF1"]

    res_dict = {}
    ## find metrics file
    for root, dirs, files in os.walk(result_dir):
        for file in files:
            if "MLP_metrics.txt" in file:
                filepath = os.path.join(root, file)
                sample_no = int(filepath.split('/')[-4])
                cell_no = int(filepath.split('/')[-2])
                if sample_no not in res_dict:
                    res_dict[sample_no] = {}
                if cell_no not in res_dict[sample_no]:
                    res_dict[sample_no][cell_no] = {}
                res_dict[sample_no][cell_no] = {}
                with open(filepath, 'r') as fopen:
                    for line in fopen:
                        metric, value = line.strip().split(":")
                        if metric in metrics:
                            res_dict[sample_no][cell_no][metric] = float(value)

    ## referred to: https://stackoverflow.com/questions/13575090/construct-pandas-dataframe-from-items-in-nested-dictionary
    df = pd.DataFrame.from_records(
        [
            (level1, level2, level3, leaf)
            for level1, level2_dict in res_dict.items()
            for level2, level3_dict in level2_dict.items()
            for level3, leaf in level3_dict.items()
        ],
        columns=['sample_seed', 'sampled_no', 'metrics', 'value']
    )

    #fig = plt.figure(figsize=(10, 6), dpi=300)
    #sns.lineplot(data=df, x="sampled_no", y="value", hue="metrics", estimator="mean", ci=95)
    #plt.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
    #if add_type == "cells":
    #    plt.xlabel("No. of cells")
    #if add_type == "inds":
    #    plt.xlabel("No. of individuals")
    #plt.ylabel("Performance")
    #plt.tight_layout()
    #plt.savefig(result_dir+os.sep+prefix+'.png')

    ## line plot on performance
    fig = plt.figure(figsize=(10, 6), dpi=300)
    clrs = sns.color_palette("tab10", len(metrics))
    for idx, metric in enumerate(metrics):
        metric_df = df[df['metrics'] == metric]
        sampled_numbers = list(set(metric_df['sampled_no']))
        sampled_numbers.sort()
        mean_list, low_qlist, high_qlist = list(), list(), list() 
        for sampled_number in sampled_numbers:
            sampled_df = metric_df[metric_df['sampled_no'] == sampled_number]
            sampled_mean = sampled_df['value'].mean()
            sampled_quantiles = mquantiles(sampled_df['value'], prob=[0.025, 0.975], 
                    alphap=1, betap=1)  ##default in R
            mean_list.append(sampled_mean)
            low_qlist.append(sampled_quantiles[0])
            high_qlist.append(sampled_quantiles[1])

        plt.plot(sampled_numbers, mean_list, label=metric, c=clrs[idx])
        plt.fill_between(sampled_numbers, 
                low_qlist, high_qlist,
                alpha=0.3, facecolor=clrs[idx])
    plt.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
    if add_type == "cells":
        plt.xlabel("No. of cells")
    if add_type == "inds":
        plt.xlabel("No. of individuals")
    plt.ylabel("Performance")
    plt.tight_layout()
    plt.savefig(result_dir+os.sep+prefix+'.png')


if __name__ == '__main__':
    ##### ================ For saturation per individual
    '''
    ### === FC sub-celltypes saturation
    res_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_saturation_collections/mousebrain_FC_sub"
    sep_index = 6
    columns = [8062, 8000, 14708, 9024, 9040, 13885]

    ### === FC major cell types saturation
    res_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_saturation_collections/mousebrain_FC_major"
    sep_index = 5
    columns = [8062, 8000, 14708, 9024, 9040, 13885]

    ### === PBMC Kang batch1 individual effect
    res_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_saturation_collections/PBMC_Kang_individualeffect"
    sep_index = 6
    columns = [1831, 1883, 1714, 1537, 1450, 2148, 1503]

    ### === PBMC Kang batch2 batch effect
    res_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_saturation_collections/PBMC_Kang_batcheffect"
    sep_index = 9
    columns = [1027, 3138, 2180, 552, 645, 2250, 2451, 2376]

    ### === PBMC Kang batch2 condition effect
    res_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_saturation_collections/PBMC_Kang_conditioneffect"
    sep_index = 9
    columns = [1363, 2703, 1998, 797, 641, 1755, 2275, 2914]

    ### === PBMC Kang batch2 combined effect
    res_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_saturation_collections/PBMC_Kang_combinedeffect"
    sep_index = 9
    columns = [2390, 5841, 4178, 1349, 1286, 4005, 4726, 5290]
    '''

    ### === mouse brain cross datasets
    #res_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_saturation_collections/mousebrain_FC_cross"
    #sep_index = 8
    ##columns = [1278, 1055, 1506, 1751, 3784, 2512, 17906, 33230, 14808]
    #columns = None

    #prefix = os.path.basename(res_dir)
    #sub_dirs = next(os.walk(res_dir))[1]
    #result_dirs = [res_dir+os.sep+x for x in sub_dirs]
    #plot_saturation_line(result_dirs, plot_dir=res_dir, sep_index=sep_index, 
    #        prefix=prefix, columns=columns)

    #### ===================== For saturation on individuals
    result_prefix = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_saturation_collections"
    result_folders = os.listdir(result_prefix)
    for fd in result_folders:
        res_dir = result_prefix+os.sep+fd
        plot_downsample_saturation_line(res_dir, prefix=fd, add_type="inds")

    ##### ====================  For saturation on downsampling
    result_prefix = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_saturation_downsample_collections"
    result_folders = os.listdir(result_prefix)
    for fd in result_folders:
        res_dir = result_prefix+os.sep+fd
        plot_downsample_saturation_line(res_dir, prefix=fd)

