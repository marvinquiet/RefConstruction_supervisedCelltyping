#!/bin/bash
#SBATCH --job-name=Mousebrain
#SBATCH --output=celltyping_datasets.out
#SBATCH --gres=gpu:1
source ~/.bashrc
conda activate /home/wma36/gpu/sc_identifier/envs/desc_tfv1

cd /homelocal/wma36/sc_identifier

PIPELINE_DIR="/home/wma36/gpu/sc_identifier/pipelines"

#declare -a celltyping_python=("MLP" "GEDFN" "ItClust" "SVM_RBF" "SVM_linear" "RF")
declare -a celltyping_python=("SVM_linear" "RF")
declare -a celltyping_R=("scmap" "CHETAH")

## === Mousebrain FC ind1 -> ind2
task="mousebrain_FC"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

## === Mousebrain FC ind1 -> ind2, using 81 sub-celltypes
task="mousebrain_FC_sub"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

## === Mousebrain FC ind1 -> ind2, using 81 sub-celltypes, purify by distance
task="mousebrain_FC_sub_purify_dist"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"


## === Mousebrain FC ind1 -> ind2, using 81 sub-celltypes, purify by SVM
task="mousebrain_FC_sub_purify_SVM"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

## Mousebrain FC exclude P60FCCx3cr1Rep1 -> P60FCCx3cr1Rep1
task="mousebrain_FC_multiinds"
train="ALL"
test="P60FCCx3cr1Rep1"

#python -m pipelines.pipeline $task -m SVM_RBF --train $train --test $test

## Mousebrain FC exclude P60FCCx3cr1Rep1 -> P60FCCx3cr1Rep1, sub-cell types
task="mousebrain_FC_multiinds_sub"
train="ALL"
test="P60FCCx3cr1Rep1"

#python -m pipelines.pipeline $task -m SVM_RBF --train $train --test $test
#Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task scmap NA NA 1000 $train $test  ## cannot run scmap result

## Mousebrain FC exclude P60FCCx3cr1Rep1 -> P60FCCx3cr1Rep1, downsample to cell numbers of P60FCRep1
task="mousebrain_FC_multiinds_sample"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

## Mousebrain FC exclude P60FCCx3cr1Rep1 -> P60FCCx3cr1Rep1, sub-cell types, downsample
task="mousebrain_FC_multiinds_sub_sample"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"


## === Mousebrain HC ind1 -> ind2
task="mousebrain_HC"
train="P60HippoRep1"
test="P60HippoRep2"

## === Mousebrain HC ind1 -> ind2, using 103 sub-celltypes
task="mousebrain_HC_sub"
train="P60HippoRep1"
test="P60HippoRep2"


## === Mousebrain FC -> HC
task="mousebrain_region"
train="FC"
test="HC"

#python -m pipelines.pipeline $task -m SVM_RBF --train $train --test $test
#Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task scmap NA NA 1000 $train $test  ## cannot run scmap result

## === Mousebrain FC -> HC, using sub-cell types
task="mousebrain_region_sub"
train="FC"
test="HC"


## === Mousebrain FC Adult -> P21
task="mousebrain_FC_stage"
train="Adult"
test="P21"


## === Mousebrain FC Adult -> P21, using sub-cell types
task="mousebrain_FC_stage_sub"
train="Adult"
test="P21"

#python -m pipelines.pipeline $task -m SVM_RBF --train $train --test $test

## === Mousebrain FC dataset1 -> dataset2
task="mousebrain_FC_datasets"
train="FC"
test="Adult"

#python -m pipelines.pipeline $task -m SVM_RBF --train $train --test $test
#python -m pipelines.pipeline $task -m RF --train $train --test $test
#Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task scmap NA NA 1000 $train $test  ## cannot run scmap result
Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task CHETAH NA NA 1000 $train $test

## === Mousebrain FC dataset1 ind1 -> dataset2 ind1
task="mousebrain_FC_datasets_multiinds"
train="FC_P60FCRep1"
test="Adult_PFCSample1"

## === Mousebrain FC dataset1 all -> dataset2 ind1, downsample to dataset1 individual
task="mousebrain_FC_datasets_multiinds_sample"
train="FC_P60FCRep1"
test="Adult_PFCSample1"

## python pipeline
#for method in ${celltyping_python[@]}; do
#    python -m pipelines.pipeline $task -m $method --train $train --test $test
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method Seurat --train $train --test $test
#    python -m pipelines.pipeline $task -m $method --select_on test --select_method Seurat --train $train --test $test
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method FEAST --train $train --test $test
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method F-test --train $train --test $test
#    python -m pipelines.pipeline $task -m $method --select_on test --select_method FEAST --train $train --test $test
#done

## R pipeline 
#for method in ${celltyping_R[@]}; do
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method NA NA 1000 $train $test
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method train Seurat 1000 $train $test
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method test Seurat 1000 $train $test
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method train FEAST 1000 $train $test
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method train F-test 1000 $train $test
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method test FEAST 1000 $train $test
#done

#for n_features in {100..5000..100} 
#do
#    echo "Select $n_features using MLP"
#    python -m pipelines.pipeline $task -m MLP --select_on train --select_method F-test --n_features $n_features --train $train --test $test
#done

conda deactivate
