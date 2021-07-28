import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
from textwrap import wrap

def extract_covariate_results(result_dir, metrics="Acc"):
    ## iterate result directory and extract metrics
    sub_dirs = next(os.walk(result_dir))[1]
    result_dirs = [result_dir+os.sep+sub_dir for sub_dir in sub_dirs]

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

def plot_metrics_together(metrics_dict, metrics="Acc", prefix=None):
    '''Plot all metrics together
    @metrics_dict: dictionary of metrics including individual metrics and other condition metrics

    ## Ref: https://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    '''
    plt.figure(figsize=(6, 6), dpi=300)

    def set_box_color(bp, color):
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)
    
    ## colors from https://htmlcolorcodes.com/
    colors = ['#CD6155', '#AF7AC5', '#5499C7', '#48C9B0', '#27AE60', '#F1C40F', '#E67E22']

    ticks = list(metrics_dict.keys())
    for i in range(len(ticks)):
        tick = ticks[i]
        box = plt.boxplot(metrics_dict[tick], positions=[i], widths=0.8)
        set_box_color(box, colors[i])
    plt.xticks(range(len(ticks)), ticks, rotation=45)
    plt.xlim(-0.5, len(ticks)-0.5)

    plt.ylabel(metrics)
    plt.tight_layout()
    plt.savefig(prefix+'_'+metrics+'.png')


if __name__ == '__main__':
    covariate_dir = "pipelines/result_TrainingCovariate_collections/"

    #prefix = "Mousebrain_Region_and_Dataset_Effect" ## comparing individual effect, dataset effect and region effect
    prefix = "PBMC_ConditonEffects_Major"

    #ind_prefix = "Mousebrain_IndividualEffect"
    #region_prefix = "Mousebrain_RegionEffect"
    ##dataset_prefix = "Mousebrain_DatasetEffect"
    #pFC_dataset_prefix = "Mousebrain_DatasetEffect_Curated"
    #cortex_dataset_prefix = "Mousebrain_DatasetEffect2"

    ind_prefix = "PBMC_IndividualEffect_Major"
    batch_prefix = "PBMC_BatchEffect_Major"
    clinical_prefix = "PBMC_ClinicalDifference_Major"

    metrics = ["Acc", "ARI", "macroF1"]
    for metric in metrics:
        metrics_dict = {}

        ## individual effects
        print("=== Individual effect:")
        result_dir = covariate_dir+os.sep+ind_prefix
        metrics_dict["Individual\nEffect"] = extract_covariate_results(result_dir,metrics=metric)

        # ==== Mouse brain datasets
        ### region effects
        #print("=== Region effect:")
        #result_dir = covariate_dir+os.sep+region_prefix
        #metrics_dict["Region\nEffect"] = extract_covariate_results(result_dir,metrics=metric)

        ### dataset effects
        #print("=== Dataset effect:")
        #result_dir = covariate_dir+os.sep+dataset_prefix
        #metrics_dict["Dataset\nEffect"] = extract_covariate_results(result_dir,metrics=metric)

        ### pFC curated dataset effects
        #print("=== preFrontal Cortex Dataset effect:")
        #result_dir = covariate_dir+os.sep+pFC_dataset_prefix
        #metrics_dict["pFC Dataset\nEffect"] = extract_covariate_results(result_dir,metrics=metric)


        ### cortex dataset effects
        #print("=== cortex dataset effect:")
        #result_dir = covariate_dir+os.sep+cortex_dataset_prefix
        #metrics_dict["cortex Dataset\nEffect"] = extract_covariate_results(result_dir,metrics=metric)

        ## ==== Human PBMC datasets

        ## Batch effect
        print("=== Batch effect:")
        result_dir = covariate_dir+os.sep+batch_prefix
        metrics_dict["Batch\nEffect"] = extract_covariate_results(result_dir,metrics=metric)

        ## Clinical Difference
        print("=== Clinical difference:")
        result_dir = covariate_dir+os.sep+clinical_prefix
        metrics_dict["Clinical\nDifference"] = extract_covariate_results(result_dir,metrics=metric)

        plot_metrics_together(metrics_dict, metrics=metric, prefix=prefix)
