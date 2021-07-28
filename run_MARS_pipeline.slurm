#!/bin/bash
#SBATCH --job-name=MARS
#SBATCH --output=celltyping_MARS.out
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=32

source ~/.bashrc

PROJECT_PATH=""
ENV_PATH=""

conda activate $ENV_PATH/mars  ## mars environment

cd $PROJECT_PATH

## 1. PBMC_batch1_indS1_to_indS5
task="PBMC_batch1_ind"
train="1154"
test="1085"

## 2. PBMC_batch1_A_to_B
task="PBMC_batch1_ABC"
train="A"
test="B"


## 3. PBMC_downsampled_batch1_A_to_indS5
task="PBMC_batch1_batchtoind"
train="A"
test="1085"

## 4. PBMC_batch2_control_to_stimulated
task="PBMC_batch2"
train="control"
test="stimulated"


## 5. PBMC_protocols_Smart-seq2_to_CEL-Seq2
task="PBMC_protocols_pbmc1"
train="Smart-seq2"
test="CEL-Seq2"


## 6. PBMC_protocols_CEL-Seq2_to_10X
task="PBMC_protocols_pbmc1"
train="CEL-Seq2"
test="10x-v2"


## 7. PBMC_protocols_10X_to_CEL-Seq2
task="PBMC_protocols_pbmc1"
train="10x-v2"
test="CEL-Seq2"

## 8. PBMC_protocols_Smart-Seq2_PBMC1_to_PBMC2
task="PBMC_protocols_batch_smart"
train="pbmc1"
test="pbmc2"


## 9. Seg_to_Muraro
task="pancreas"
train="seg"
test="muraro"


## 10. Muraro_to_Seg
task="pancreas"
train="muraro"
test="seg"


## 11. Seg_Xin_to_Muraro
task="pancreas_multi_to_multi"
train="seg_xin"
test="muraro"

## 12. Seg_Healthy_to_T2D
task="pancreas_seg_cond"
train="Healthy"
test="T2D"

## 13. Seg_Healthy_to_Muraro
task="pancreas_custom"
train="Healthy"
test="muraro"


## 14. Seg_6Healthy3T2D_to_Seg_1T2D
task="pancreas_seg_mix"
train="Healthy"
test="T2D"


## 15. Seg_5Healthy4T2D_to_Seg_1Healthy
task="pancreas_seg_mix"
train="T2D"
test="Healthy"


## 16. Wholebrain_FC_ind1_to_ind2_major_celltypes
task="mousebrain_FC"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

## 17. Wholebrain_FC_ind1_to_ind2_sub_celltypes
task="mousebrain_FC_sub"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"


## 18. Wholebrain_HC_ind1_to_ind2_major_celltypes
task="mousebrain_HC"
train="P60HippoRep1"
test="P60HippoRep2"


## 19. Wholebrain_HC_ind1_to_ind2_sub_celltypes
task="mousebrain_HC_sub"
train="P60HippoRep1"
test="P60HippoRep2"

## 20. Wholebrain_FC_to_HC_major_celltypes, cannot run no-feature because CUDA out of memory
task="mousebrain_region"
train="FC"
test="HC"


## 21. pFC_Adult_to_P21_major_celltypes
task="mousebrain_FC_stage"
train="Adult"
test="P21"


## 22. pFC_Adult_to_P21_sub_celltypes
task="mousebrain_FC_stage_sub"
train="Adult"
test="P21"


## 23. Wholebrain_FC_adult_to_pFC_adult_major_celltypes, cannot run no-feature because CUDA out of memory
task="mousebrain_FC_datasets"
train="FC"
test="Adult"


## 24. Wholebrain_FC_ind1_to_pFC_ind1
task="mousebrain_FC_datasets_multiinds"
train="FC_P60FCRep1"
test="Adult_PFCSample1"


## 25. Wholebrain_FC_exclude_ind2_to_ind2_major_celltypes
task="mousebrain_FC_multiinds"
train="ALL"
test="P60FCCx3cr1Rep1"


## 26. Wholebrain_FC_exclude_ind2_to_ind2_sub_celltypes
task="mousebrain_FC_multiinds_sub"
train="ALL"
test="P60FCCx3cr1Rep1"


## 27. Wholebrain_FC_downsampled_exclude_ind2_to_ind2_major_celltypes
task="mousebrain_FC_multiinds_sample"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

## 28. Wholebrain_FC_downsample_exclude_ind2_to_ind2_sub_celltypes
task="mousebrain_FC_multiinds_sub_sample"
train="P60FCRep1"
test="P60FCCx3cr1Rep1"

## 29. Wholebrain_FC_downsample_adult_to_pFC_ind1
task="mousebrain_FC_datasets_multiinds_sample"
train="FC_P60FCRep1"
test="Adult_PFCSample1"

### ============================== separator for revision

## 30. PBMC_Zheng_FACS cross validation
task="PBMC_Zheng_FACS"
train="0.8"
test="0.2"

## 31. PBMC_Zheng_FACS cross validation on major cell types
task="PBMC_Zheng_FACS_curated"
train="0.8"
test="0.2"

## 32. PBMC cross dataset, Kang -> Zheng (0.2 target)
task="PBMC_cross"
train="Kang"
test="Zheng"

## 33. PBMC cross dataset, Ding -> Zheng (0.2 target)
task="PBMC_cross"
train="Ding"
test="Zheng"

## 34. Allen brain Smart-Seqv4 cross-validation, 80% train -> 20% test
task="allenbrain_ss"
train="0.8"
test="0.2"

## 35. Allen brain 10X ind1 -> ind2
task="allenbrain_10x"
train="372312"
test="371230"

## 36. PBMC cross dataset, Kang batch1 -> Zheng (0.2 target)
task="PBMC_cross"
train="Kang_batch1"
test="Zheng"

## 36. PBMC cross dataset, Ding droplet -> Zheng (0.2 target)
task="PBMC_cross"
train="Ding_droplet"
test="Zheng"

## 37. mouse FC dataset -> allen brain SSv4 (0.2 target)
task="allenbrain_cross"
train="FC"
test="0.2"

## 38. mouse pFC dataset -> allen brain SSv4 (0.2 target)
task="allenbrain_cross"
train="pFC"
test="0.2"

python -m test_MARS.test_MARS $task --train $train --test $test
python -m test_MARS.test_MARS $task --select_on train --select_method Seurat --train $train --test $test
python -m test_MARS.test_MARS $task --select_on test --select_method Seurat --train $train --test $test
python -m test_MARS.test_MARS $task --select_on train --select_method FEAST --train $train --test $test
python -m test_MARS.test_MARS $task --select_on train --select_method F-test --train $train --test $test
python -m test_MARS.test_MARS $task --select_on test --select_method FEAST --train $train --test $test


