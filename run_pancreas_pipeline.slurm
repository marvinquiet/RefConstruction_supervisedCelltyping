#!/bin/bash
#SBATCH --job-name=Pancreas
#SBATCH --output=celltyping_Pancreas.out
#SBATCH --gres=gpu:1
source ~/.bashrc

PROJECT_PATH=""
ENV_PATH=""

conda activate $ENV_PATH/celltyping  

cd $PROJECT_PATH

PIPELINE_DIR=$PROJECT_PATH/pipelines

declare -a celltyping_python=("MLP" "GEDFN" "ItClust" "SVM_RBF" "SVM_linear" "RF")
declare -a celltyping_R=("scmap" "CHETAH")

## === Seg -> Muraro
task="pancreas"
train="seg"
test="muraro"


## === Muraro -> Seg
task="pancreas"
train="muraro"
test="seg"

## === Seg Healthy -> T2D
task="pancreas_seg_cond"
train="Healthy"
test="T2D"


## === Seg Healthy -> Muraro
task="pancreas_custom"
train="Healthy"
test="muraro"


## Seg Healthy+T2D -> T2D
task="pancreas_seg_mix"
train="Healthy"
test="T2D"


## Seg Healthy+T2D -> Healthy
task="pancreas_seg_mix"
train="T2D"
test="Healthy"


## Seg+Xin -> Muraro
task="pancreas_multi_to_multi"
train="seg_xin"
test="muraro"


## python pipeline
#for method in ${celltyping_python[@]}; do
#    python -m pipelines.pipeline $task -m $method --train $train --test $test
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method Seurat --train $train --test $test
#    python -m pipelines.pipeline $task -m $method --select_on test --select_method Seurat --train $train --test $test
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method FEAST --train $train --test $test
#    python -m pipelines.pipeline $task -m $method --select_on test --select_method FEAST --train $train --test $test
#done

## R pipeline
#for method in ${celltyping_R[@]}; do
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method NA NA 1000 $train $test
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method train Seurat 1000 $train $test
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method test Seurat 1000 $train $test
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method train FEAST 1000 $train $test
#    Rscript $PIPELINE_DIR/Rcode/scmap_CHETAH_pipeline.R $task $method test FEAST 1000 $train $test
#done


conda deactivate
