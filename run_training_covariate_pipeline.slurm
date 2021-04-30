#!/bin/bash
#SBATCH --job-name=TrainingCovariate
#SBATCH --output=training_covariate_MouseBothCortex.out
#SBATCH --gres=gpu:1
source ~/.bashrc
conda activate /home/wma36/gpu/sc_identifier/envs/desc_tfv1

cd /homelocal/wma36/sc_identifier

PIPELINE_DIR="/home/wma36/gpu/sc_identifier/pipelines"

## fix method and feature selection
method="MLP"
FS="F-test"

### === cross region, use mouse brain 
## set individual 2 as fixed target
test="P60FCCx3cr1Rep1"
declare -a HC_train_list=("P60HippoRep1" "P60HippoRep2" "P60HippoRep3" "P60HippoRep4" "P60HippoRep5" "P60HippoRep6")

task="mousebrain_region"
## HC ind1-6 -> FC ind2
#for train in ${HC_train_list[@]}; do
#        python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train "HC_${train}" --test "FC_${test}"
#done

## use all HC individuals to run
#python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train "HC" --test "FC_${test}"

### === cross dataset, use mouse brain pFC
## same predictor
test="P60FCCx3cr1Rep1"
declare -a pFC_train_list=("2" "3" "5" "7" "11" "12") ## saline samples, take as adult controls

task="mousebrain_FC_datasets_multiinds"
## pFC inds -> FC ind2
#for train in ${pFC_train_list[@]}; do
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train "Adult_PFCSample${train}" --test "FC_${test}"
#done

## use all pFC inds to run
#python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train "Adult-Saline" --test "FC_${test}"

### === cross dataset, use mouse cortex protocol
## same predictor
test="P60FCCx3cr1Rep1"
declare -a proc_train_list=("+Cortex1+DroNc-seq" "+Cortex2+DroNc-seq")

task="mousecortex_protocols_mouseFC_ind"
#for train in ${proc_train_list[@]}; do
#    python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train "MouseProtocol${train}" --test "mouseFC_${test}"
#done
## use both cortex1+cortex2 to predict
#python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train "MouseProtocol+Both+DroNc-seq" --test "mouseFC_${test}"


## === cross dataset, using all mouse brain datasets -> FC ind1
task="mousebrain_combined_mouseFC"
test="P60FCCx3cr1Rep1"
#python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train "MouseProtocol+Both+DroNc-seq_pFC+Adult+Saline" --test "FC_${test}"

### === cross condition, use Human PBMC
# set individual S5 as fixed target
test="1085"
declare -a PBMC_batch1_train_list=("1043" "1079" "1154" "1249" "1493" "1511" "1598") ## already have this result
task="PBMC_batch1_ind_major"
for train in ${PBMC_batch1_train_list[@]}; do
    python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train $train --test $test
done


declare -a PBMC_batch2_train_list=("101" "1015" "1016" "1039" "107" "1244" "1256" "1488")
#task="PBMC_batch1_batch2_ind"
task="PBMC_batch1_batch2_ind_major"
for train in ${PBMC_batch2_train_list[@]}; do
    python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train "batch2ctrl_${train}" --test "batch1_${test}"
    python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train "batch2stim_${train}" --test "batch1_${test}"
done

python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train batch2ctrl --test "batch1_${test}"
python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train batch2stim --test "batch1_${test}"
python -m pipelines.pipeline $task -m $method --select_on train --select_method $FS --train batch2 --test "batch1_${test}"
