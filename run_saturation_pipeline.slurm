#!/bin/bash
#SBATCH --job-name=saturation
#SBATCH --output=celltyping_datasets_saturation.out
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=32

source ~/.bashrc
PROJECT_PATH=""
ENV_PATH=""

conda activate $ENV_PATH/celltyping  

cd $PROJECT_PATH

declare -a sample_seeds=(5866 9558 3578 8268 2281 4617 2289 1553 4104 8725 9861 2407 5081 1618 1208 5409 7735 9171 1649 5796 7113 5180 3350 9052 7815 7253 8541 4267 1020 8989 230 1528 6534 18 8086 5458 3996 5328 1031 3130) ## sampled by random.sample(range(0, 10000), k=10) with random.seed(0)
#declare -a sample_seeds=(0 6311 6890 663 4242 8376 7961 6634 4969 7808)

# =========== Add by individual

## === Mouse brain FC saturation test, test for major/sub-cell types
#task="mousebrain_FC_sub"
#test="P60FCCx3cr1Rep1"
#declare -a inds=("P60FCRep2" "P60FCRep3" "P60FCRep4" "P60FCRep6" "P60FCAldh1l1Rep1")
#train="P60FCRep1"
#for ind in ${inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --sample_seed $seed
#done
#
#task="mousebrain_FC"
#test="P60FCCx3cr1Rep1"
#declare -a inds=("P60FCRep2" "P60FCRep3" "P60FCRep4" "P60FCRep6" "P60FCAldh1l1Rep1")
#train="P60FCRep1"
#for ind in ${inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --sample_seed $seed
#done
#
### === Mouse brain cross datasets -> mouse FC target
#task="mousebrain_crossdataset_inds"
#test="mFC_P60FCCx3cr1Rep1"
#declare -a pFC_inds=("PFCSample3" "PFCSample5" "PFCSample7" "PFCSample11" "PFCSample12")
#declare -a allen_inds=("371230" "372312" "372314")
#train="pFC_PFCSample2"
#for ind in ${pFC_inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#train+="+allenbrain"
#for ind in ${allen_inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --sample_seed $seed
#done


### === PBMC Kang combine batch1 individuals -> batch1 target
#task="PBMC_batch1_ind"
#test="1085"
#declare -a inds=("1043" "1079" "1249" "1493" "1511" "1598")
#train="1154"
#for ind in ${inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --sample_seed $seed
#done
#
### === PBMC Kang combine batch2 indivdiuals from same condition -> batch1 target
#task="PBMC_batch1_batch2_ind"
#test="batch1_1085"
#declare -a inds=("1015" "1016" "1039" "107" "1244" "1256" "1488")
#train="batch2control_101"
#for ind in ${inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --sample_seed $seed
#done
#
#
### === PBMC Kang combine batch2 individuals from different condition -> batch1 target
#task="PBMC_batch1_batch2_ind"
#test="batch1_1085"
#declare -a inds=("1015" "1016" "1039" "107" "1244" "1256" "1488")
#train="batch2stimulated_101"
#for ind in ${inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --sample_seed $seed
#done
#
#
### === PBMC Kang combine batch2 from two conditions -> batch1 target
#task="PBMC_batch1_batch2_ind"
#test="batch1_1085"
#declare -a inds=("1015" "1016" "1039" "107" "1244" "1256" "1488")
#train="batch2_101"
#for ind in ${inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --sample_seed $seed
#done


## ============ Combine all samples and do downsample
## === Mouse brain FC saturation test, test for major/sub-cell types
task="mousebrain_FC_sub"
test="P60FCCx3cr1Rep1"
declare -a inds=("P60FCRep2" "P60FCRep3" "P60FCRep4" "P60FCRep6" "P60FCAldh1l1Rep1")
train="P60FCRep1"
for ind in ${inds[@]}; do
    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
done

for seed in ${sample_seeds[@]}; do
    python -m pipelines.saturation_pipeline $task  --train $train --test $test --downsample --sample_seed $seed
done

task="mousebrain_FC"
test="P60FCCx3cr1Rep1"
declare -a inds=("P60FCRep2" "P60FCRep3" "P60FCRep4" "P60FCRep6" "P60FCAldh1l1Rep1")
train="P60FCRep1"
for ind in ${inds[@]}; do
    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
done

for seed in ${sample_seeds[@]}; do
    python -m pipelines.saturation_pipeline $task  --train $train --test $test --downsample --sample_seed $seed
done

## === Mouse brain cross datasets -> mouse FC target
task="mousebrain_crossdataset_inds"
test="mFC_P60FCCx3cr1Rep1"
declare -a pFC_inds=("PFCSample3" "PFCSample5" "PFCSample7" "PFCSample11" "PFCSample12")
declare -a allen_inds=("371230" "372312" "372314")
train="pFC_PFCSample2"
for ind in ${pFC_inds[@]}; do
    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
done
train+="+allenbrain"
for ind in ${allen_inds[@]}; do
    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
done

for seed in ${sample_seeds[@]}; do
    python -m pipelines.saturation_pipeline $task  --train $train --test $test --downsample --sample_seed $seed
done

## === PBMC Kang combine batch1 individuals -> batch1 target
#task="PBMC_batch1_ind"
#test="1085"
#declare -a inds=("1043" "1079" "1249" "1493" "1511" "1598")
#train="1154"
#for ind in ${inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --downsample --sample_seed $seed --downsample_size 1000
#done
#
### === PBMC Kang combine batch2 indivdiuals from same condition -> batch1 target
#task="PBMC_batch1_batch2_ind"
#test="batch1_1085"
#declare -a inds=("1015" "1016" "1039" "107" "1244" "1256" "1488")
#train="batch2control_101"
#for ind in ${inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --downsample --sample_seed $seed --downsample_size 1000
#done
#
### === PBMC Kang combine batch2 individuals from different condition -> batch1 target
#task="PBMC_batch1_batch2_ind"
#test="batch1_1085"
#declare -a inds=("1015" "1016" "1039" "107" "1244" "1256" "1488")
#train="batch2stimulated_101"
#for ind in ${inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --downsample --sample_seed $seed --downsample_size 1000
#done
#
#
### === PBMC Kang combine batch2 from two conditions -> batch1 target
#task="PBMC_batch1_batch2_ind"
#test="batch1_1085"
#declare -a inds=("1015" "1016" "1039" "107" "1244" "1256" "1488")
#train="batch2_101"
#for ind in ${inds[@]}; do
#    train=$(paste <(echo "$train") <(echo "$ind") -d '_')
#done
#
#for seed in ${sample_seeds[@]}; do
#    python -m pipelines.saturation_pipeline $task  --train $train --test $test --downsample --sample_seed $seed --downsample_size 1000
#done

conda deactivate
