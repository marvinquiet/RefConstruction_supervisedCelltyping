#!/bin/bash
#SBATCH --job-name=twostage
#SBATCH --output=celltyping_twostage_mousebrain.out
#SBATCH --gres=gpu:1
source ~/.bashrc
conda activate /home/wma36/gpu/sc_identifier/envs/desc_tfv1

cd /homelocal/wma36/sc_identifier

PIPELINE_DIR="/home/wma36/gpu/sc_identifier/pipelines"

declare -a celltyping_python=("MLP" "GEDFN" "ItClust" "SVM_RBF" "SVM_linear" "RF")

## === Mousebrain FC ind1 -> ind2, two-stage
task="mousebrain_FC"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

## python pipeline
for method in ${celltyping_python[@]}; do
    python -m pipelines.twostage_pipeline $task -m $method --select_on train --select_method F-test --train $train --test $test
done


## === PBMC batch1 A -> B, two-stage
task="PBMC_batch1"
train="A"
test="B"

## === PBMC batch1 A -> B, two-stage
task="PBMC_batch2"
train="control"
test="stimulated"

## python pipeline
for method in ${celltyping_python[@]}; do
    #python -m pipelines.twostage_pipeline $task -m $method --train $train --test $test
    #python -m pipelines.twostage_pipeline $task -m $method --select_on train --select_method Seurat --train $train --test $test
    #python -m pipelines.twostage_pipeline $task -m $method --select_on test --select_method Seurat --train $train --test $test
    #python -m pipelines.twostage_pipeline $task -m $method --select_on train --select_method FEAST --train $train --test $test
    #python -m pipelines.twostage_pipeline $task -m $method --select_on train --select_method F-test --train $train --test $test
    #python -m pipelines.twostage_pipeline $task -m $method --select_on test --select_method FEAST --train $train --test $test
done

conda deactivate
