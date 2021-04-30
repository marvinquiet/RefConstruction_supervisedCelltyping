#!/bin/bash
#SBATCH --job-name=PBMC
#SBATCH --output=celltyping_PBMC.out
#SBATCH --gres=gpu:1
source ~/.bashrc
conda activate /home/wma36/gpu/sc_identifier/envs/desc_tfv1

cd /homelocal/wma36/sc_identifier

PIPELINE_DIR="/home/wma36/gpu/sc_identifier/pipelines"

declare -a celltyping_python=("MLP" "GEDFN" "ItClust" "SVM_RBF" "SVM_linear" "RF")
#declare -a celltyping_python=("MLP_GO" "MLP_CP")
declare -a celltyping_R=("scmap" "CHETAH")

## ==== PBMC batch1 S1 -> S5
task="PBMC_batch1_ind"
train="1154"
test="1085"

## python pipeline
# **jiaying**: using CUDA_VISIBLE_DEVICES=-1 and 
# directly `sh thisfile.slurm` instead of `slurm thisfile.slurm`
#CUDA_VISIBLE_DEVICES=-1 python -m pipelines.pipeline $task -m MLP --select_on train --select_method F-test --train $train --test $test
#CUDA_VISIBLE_DEVICES=-1 python -m pipelines.pipeline $task -m GEDFN --select_on train --select_method F-test --train $train --test $test
#CUDA_VISIBLE_DEVICES=-1 python -m pipelines.pipeline $task -m ItClust --select_on train --select_method F-test --train $train --test $test

## ==== PBMC batch1 A -> B
task="PBMC_batch1_ABC"
train="A"
test="B"


## === PBMC batch2 control -> stimulated
task="PBMC_batch2"
train="control"
test="stimulated"

## python pipeline
# **jiaying**: using CUDA_VISIBLE_DEVICES=-1 and 
# directly `sh thisfile.slurm` instead of `slurm thisfile.slurm`
CUDA_VISIBLE_DEVICES=-1 python -m pipelines.pipeline $task -m MLP --select_on train --select_method F-test --train $train --test $test
CUDA_VISIBLE_DEVICES=-1 python -m pipelines.pipeline $task -m GEDFN --select_on train --select_method F-test --train $train --test $test
CUDA_VISIBLE_DEVICES=-1 python -m pipelines.pipeline $task -m ItClust --select_on train --select_method F-test --train $train --test $test


## ==== PBMC batch1 A -> B purify A using distance
task="PBMC_batch1_ABC_purify_dist"
train="A"
test="B"

## ==== PBMC batch1 A -> B purify A using SVM
task="PBMC_batch1_ABC_purify_SVM"
train="A"
test="B"

## === PBMC batch2 control -> stimulated purify control using distance
task="PBMC_batch2_purify_dist"
train="control"
test="stimulated"

## === PBMC batch2 control -> stimulated purify control using SVM
task="PBMC_batch2_purify_SVM"
train="control"
test="stimulated"

## === PBMC batch1 S1 -> S5 purify control using distance
task="PBMC_batch1_ind_purify_dist"
train="1154"
test="1085"

## === PBMC batch1 S1 -> S5 purify control using SVM
task="PBMC_batch1_ind_purify_SVM"
train="1154"
test="1085"

## === PBMC batch1 A downsampled  -> S5
task="PBMC_batch1_batchtoind"
train="A"
test="1085"

## === PBMC protocols plate-based, pbmc1, Smart-seq2 -> CEL-Seq2
task="PBMC_protocols_pbmc1"
train="Smart-seq2"
test="CEL-Seq2"

## === PBMC protocols pbmc1, 10X -> CEL-seq2
task="PBMC_protocols_pbmc1"
train="10x-v2"
test="CEL-Seq2"

## === PBMC protocols pbmc1, 10X -> CEL-seq2
task="PBMC_protocols_pbmc1"
train="CEL-Seq2"
test="10x-v2"

## === PBMC protocols plate-based, Smart-seq2,  pbmc1-> pbmc2
task="PBMC_protocols_batch_smart"
train="pbmc1"
test="pbmc2"

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


### ======= for different feature numbers
#for n_features in {100..5000..100} 
#do
#    echo "Select $n_features using MLP"
#    python -m pipelines.pipeline $task -m MLP --select_on train --select_method FEAST_F-test --n_features $n_features --train $train --test $test
#done



conda deactivate
