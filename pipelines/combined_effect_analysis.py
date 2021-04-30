import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
from textwrap import wrap

def extract_combinedeffect_results(result_dirs, metrics="Acc"):
    metrics_list = []
    for result_dir in result_dirs:
        metrics_path = result_dir+os.sep+'F-test_1000_on_train'+os.sep+'MLP_metrics.txt'
        with open(metrics_path, 'r') as fopen:
            for line in fopen:
                metric, value = line.split(":")
                if metric == metrics:
                    metrics_list.append(float(value))

    print(metrics)
    print("Mean:", np.mean(metrics_list))
    print("Variance:", np.var(metrics_list))

    return metrics_list

def plot_three_together(individual_metrics, downsample_metrics, all_metrics, 
        metrics="Acc", prefix=None):
    '''Plot all three metrics together using all_metrics as a red line

    ## Ref: https://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    '''
    plt.figure(figsize=(6, 6), dpi=300)

    def set_box_color(bp, color):
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)

    ticks = ["Individual", "Combined\n Individual\n Downsampled"]
    box_ind = plt.boxplot(individual_metrics, positions=[0], widths=0.8)
    box_downsample = plt.boxplot(downsample_metrics, positions=[1], widths=0.8)
    set_box_color(box_ind, '#D7191C') # colors are from http://colorbrewer2.org/
    set_box_color(box_downsample, '#2C7BB6')

    if len(all_metrics) == 1:
        ## draw line of combining all metrics together
        plt.axhline(y=all_metrics, color="k", linestyle='-')

        plt.xticks([0, 1], ticks)
        plt.xlim(-0.5, 1.5)
    else:
        ticks.append("Combined\n Individual") 
        box_all = plt.boxplot(all_metrics, positions=[2], widths=0.8)
        set_box_color(box_all, "#000000")

        plt.xticks([0, 1, 2], ticks)
        plt.xlim(-0.5, 2.5)

    plt.ylabel(metrics)
    plt.tight_layout()
    plt.savefig(prefix+'_'+metrics+'.png')


if __name__ == '__main__':
    combined_effect_dir = "result_CombinedEffect_collections/"

    #prefix = "Mouse_brain_interdataset_CombinedEffect"  ## Human_PBMC_CombinedEffect, Mouse_brain_major_CombinedEffect
    ## for Mouse_brain_interdataset_CombinedEffect, the ALL -> one metrics is also a boxplot

    ## all prefix
    #all_prefix = "result_mousebrain_FC_datasets_ALL_multiindis_pFC_adult_inds"
    ### downsample prefix
    #downsample_prefix = "result_mousebrain_FC_datasets_sample_Avg_multiindis_pFC_adult_inds"
    ### individual prefix
    #ind_prefix = "result_mousebrain_FC_datasets_multiindis_pFC_adult_inds"

    #prefix = "Human_PBMC_CombinedEffect"
    #all_prefix = "result_PBMC_batch1_multiinds_ALL_to_1085"
    #downsample_prefix = "result_PBMC_batch1_multiinds_sample_ALL_to_1085"
    #ind_prefix = "result_PBMC_batch1_inds"

    #prefix = "Mouse_brain_intradataset_CombinedEffect_major"
    #all_prefix = "result_mousebrain_FC_multiinds_ALL_to_P60FCCx3cr1Rep1"
    #downsample_prefix = "result_mousebrain_FC_multiinds_sample_Avg_to_P60FCCx3cr1Rep1"
    #ind_prefix = "result_mousebrain_FC_inds"

    prefix = "Mouse_brain_intradataset_CombinedEffect_sub"
    all_prefix = "result_mousebrain_FC_multiinds_sub_ALL_to_P60FCCx3cr1Rep1"
    downsample_prefix = "result_mousebrain_FC_multiinds_sub_sample_Avg_to_P60FCCx3cr1Rep1"
    ind_prefix = "result_mousebrain_FC_sub_inds"

    metrics = ["Acc", "ARI", "macroF1"]
    for metric in metrics:
        print("===Downsample:")
        ## downsample metrics
        result_dir = combined_effect_dir+os.sep+downsample_prefix
        sub_dirs = next(os.walk(result_dir))[1]
        result_dirs = [result_dir+os.sep+sub_dir for sub_dir in sub_dirs]
        downsample_metrics_list = extract_combinedeffect_results(result_dirs,
                metrics=metric)

        print("===Individual:")
        ## individual metrics
        result_dir = combined_effect_dir+os.sep+ind_prefix
        sub_dirs = next(os.walk(result_dir))[1]
        result_dirs = [result_dir+os.sep+sub_dir for sub_dir in sub_dirs]
        individual_metrics_list = extract_combinedeffect_results(result_dirs, 
                metrics=metric)

        print("===All to one:")
        ## all metrics to one
        result_dir = combined_effect_dir+os.sep+all_prefix
        sub_dirs = next(os.walk(result_dir))[1]
        #result_dirs = [result_dir+os.sep+sub_dir for sub_dir in sub_dirs]
        result_dirs = [result_dir]
        all_metrics_list = extract_combinedeffect_results(result_dirs, metrics=metric)

        ## extract all metric
        #all_metrics_path = combined_effect_dir+os.sep+all_prefix+os.sep+'F-test_1000_on_train'+os.sep+'MLP_metrics.txt'
        #all_metrics = None
        #with open(all_metrics_path, 'r') as fopen:
        #    for line in fopen:
        #        met, value = line.split(":")
        #        if met == metric:
        #            all_metrics = float(value)
        #print("===All metrics:", all_metrics)

        ## plot these three information together
        plot_three_together(individual_metrics_list, downsample_metrics_list, all_metrics_list,
                metrics=metric, prefix=prefix)

