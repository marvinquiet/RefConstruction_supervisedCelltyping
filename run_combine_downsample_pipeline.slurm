#!/bin/bash
#SBATCH --job-name=Downsample
#SBATCH --output=celltyping_combine_effect_MousebrainCurated.out
#SBATCH --gres=gpu:1
source ~/.bashrc
conda activate /home/wma36/gpu/sc_identifier/envs/desc_tfv1

cd /homelocal/wma36/sc_identifier

PIPELINE_DIR="/home/wma36/gpu/sc_identifier/pipelines"

## fix method and feature selection
method="MLP"
FS="F-test"

##### =============== Mouse Brain Dataset
## set individual 2 as fixed target
test="P60FCCx3cr1Rep1"
declare -a train_list=("P60FCRep1" "P60FCRep2" "P60FCRep3" "P60FCRep4" "P60FCRep6" "P60FCAldh1l1Rep1")

## === Mousebrain FC each ind1 -> ind2
task="mousebrain_FC"

## === Mouse FC each ind1 -> ind2, with cell types curated (layers neuron -> neuron; two interneuron -> interneuron)
task="mousebrain_FC_curate"
#for train in ${train_list[@]}; do
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train $train --test $test
#done


## === Mousebrain FC ind1 -> ind2, using 81 sub-celltypes
task="mousebrain_FC_sub"

## Mousebrain FC exclude P60FCCx3cr1Rep1 -> P60FCCx3cr1Rep1
task="mousebrain_FC_multiinds"
train="ALL"

## Mousebrain FC exclude P60FCCx3cr1Rep1 -> P60FCCx3cr1Rep1, sub-cell types
task="mousebrain_FC_multiinds_sub"
train="ALL"

## Mousebrain FC exclude P60FCCx3cr1Rep1 -> P60FCCx3cr1Rep1, downsample to train cell numbers
task="mousebrain_FC_multiinds_sample"
train="Avg"

## Mousebrain FC exclude P60FCCx3cr1Rep1 -> P60FCCx3cr1Rep1, sub-cell types, downsample
task="mousebrain_FC_multiinds_sub_sample"
train="Avg"

#### ================ Mouse Brain across dataset
## TODO: wait on later
declare -a FC_train_list=("P60FCRep1" "P60FCRep2" "P60FCRep3" "P60FCRep4" "P60FCRep6" "P60FCAldh1l1Rep1" "P60FCCx3cr1Rep1") ## adult controls
declare -a pFC_test_list=("2" "3" "5" "7" "11" "12") ## saline samples, take as adult controls

## Mouse brain FC individual -> pFC Adult individual
task="mousebrain_FC_datasets_multiinds"
#for train in ${FC_train_list[@]}; do
#    for test in ${pFC_test_list[@]}; do
#        python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train "FC_${train}" --test "Adult_PFCSample${test}"
#    done
#done

## Mouse brain all FC individual -> pFC Adult individuals 
task="mousebrain_FC_datasets_multiinds"
train="FC"
#for test in ${pFC_test_list[@]}; do
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train $train --test "Adult_PFCSample${test}"
#done

task="mousebrain_FC_datasets_multiinds_sample"
train="FC"

#for i in {1..11}
#do
#    RANDOM_NUM=$(($RANDOM % 1000))
#    echo "Random downsampling seed, $RANDOM_NUM"
#    for test in ${pFC_test_list[@]}; do
#        python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train $train --test "Adult_PFCSample${test}" --sample_seed $RANDOM_NUM
#    done
#done



##### =============== Human PBMC Dataset
## set individual S5 as fixed target
test="1085"
declare -a PBMC_train_list=("1043" "1079" "1154" "1249" "1493" "1511" "1598")

## Human PBMC batch1 ind1 -> ind2
task="PBMC_batch1_ind"
#for train in ${PBMC_train_list[@]}; do
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train $train --test $test
#done

## Human PBMC batch1 exclude ind2 -> ind2
task="PBMC_batch1_multiinds"
train="ALL"


## Human PBMC batch1 exclude ind2 -> ind2, downsample
task="PBMC_batch1_multiinds_sample"
train="ALL"


#for i in {1..30}
#do
#    RANDOM_NUM=$(($RANDOM % 1000))
#    echo "Random downsampling seed, $RANDOM_NUM"
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train $train --test $test --sample_seed $RANDOM_NUM
#done


conda deactivate
