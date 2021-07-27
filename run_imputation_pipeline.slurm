#!/bin/bash
#SBATCH --job-name=imputation
#SBATCH --output=celltyping_imputation_SAVER.out
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10

source ~/.bashrc
conda activate /home/wma36/gpu/envs/celltyping  ## mars environment

cd /home/wma36/gpu/celltyping_refConstruct/

#imputation="MAGIC"
#imputation="scVI"
imputation="SAVER"

## 1. PBMC_batch1_indS1_to_indS5
task="PBMC_batch1_ind"
train="1154"
test="1085"

#python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 2. PBMC_batch1_A_to_B
task="PBMC_batch1_ABC"
train="A"
test="B"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

exit 0

## 3. PBMC_downsampled_batch1_A_to_indS5
task="PBMC_batch1_batchtoind"
train="A"
test="1085"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 4. PBMC_batch2_control_to_stimulated
task="PBMC_batch2"
train="control"
test="stimulated"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 5. PBMC_protocols_Smart-seq2_to_CEL-Seq2
task="PBMC_protocols_pbmc1"
train="Smart-seq2"
test="CEL-Seq2"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 6. PBMC_protocols_CEL-Seq2_to_10X
task="PBMC_protocols_pbmc1"
train="CEL-Seq2"
test="10x-v2"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 7. PBMC_protocols_10X_to_CEL-Seq2
task="PBMC_protocols_pbmc1"
train="10x-v2"
test="CEL-Seq2"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 8. PBMC_protocols_Smart-Seq2_PBMC1_to_PBMC2
task="PBMC_protocols_batch_smart"
train="pbmc1"
test="pbmc2"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 9. Seg_to_Muraro
task="pancreas"
train="seg"
test="muraro"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 10. Muraro_to_Seg
task="pancreas"
train="muraro"
test="seg"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 11. Seg_Xin_to_Muraro
task="pancreas_multi_to_multi"
train="seg_xin"
test="muraro"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 12. Seg_Healthy_to_T2D
task="pancreas_seg_cond"
train="Healthy"
test="T2D"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 13. Seg_Healthy_to_Muraro
task="pancreas_custom"
train="Healthy"
test="muraro"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 14. Seg_6Healthy3T2D_to_Seg_1T2D
task="pancreas_seg_mix"
train="Healthy"
test="T2D"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 15. Seg_5Healthy4T2D_to_Seg_1Healthy
task="pancreas_seg_mix"
train="T2D"
test="Healthy"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 16. Wholebrain_FC_ind1_to_ind2_major_celltypes
task="mousebrain_FC"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 17. Wholebrain_FC_ind1_to_ind2_sub_celltypes
task="mousebrain_FC_sub"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 18. Wholebrain_HC_ind1_to_ind2_major_celltypes
task="mousebrain_HC"
train="P60HippoRep1"
test="P60HippoRep2"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 19. Wholebrain_HC_ind1_to_ind2_sub_celltypes
task="mousebrain_HC_sub"
train="P60HippoRep1"
test="P60HippoRep2"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 20. Wholebrain_FC_to_HC_major_celltypes, cannot run no-feature because CUDA out of memory
task="mousebrain_region"
train="FC"
test="HC"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 21. pFC_Adult_to_P21_major_celltypes
task="mousebrain_FC_stage"
train="Adult"
test="P21"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 22. pFC_Adult_to_P21_sub_celltypes
task="mousebrain_FC_stage_sub"
train="Adult"
test="P21"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 23. Wholebrain_FC_adult_to_pFC_adult_major_celltypes, cannot run no-feature because CUDA out of memory
task="mousebrain_FC_datasets"
train="FC"
test="Adult"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 24. Wholebrain_FC_ind1_to_pFC_ind1
task="mousebrain_FC_datasets_multiinds"
train="FC_P60FCRep1"
test="Adult_PFCSample1"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 25. Wholebrain_FC_exclude_ind2_to_ind2_major_celltypes
task="mousebrain_FC_multiinds"
train="ALL"
test="P60FCCx3cr1Rep1"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 26. Wholebrain_FC_exclude_ind2_to_ind2_sub_celltypes
task="mousebrain_FC_multiinds_sub"
train="ALL"
test="P60FCCx3cr1Rep1"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 27. Wholebrain_FC_downsampled_exclude_ind2_to_ind2_major_celltypes
task="mousebrain_FC_multiinds_sample"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 28. Wholebrain_FC_downsample_exclude_ind2_to_ind2_sub_celltypes
task="mousebrain_FC_multiinds_sub_sample"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation

## 29. Wholebrain_FC_downsample_adult_to_pFC_ind1
task="mousebrain_FC_datasets_multiinds_sample"
train="FC_P60FCRep1"
test="Adult_PFCSample1"

python -m test_imputation.test_imputation $task --train $train --test $test --imputation $imputation
