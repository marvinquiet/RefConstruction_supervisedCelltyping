# Evaluation of some aspects in supervised cell type identification for single-cell RNA-seq: classifier, feature selection, and reference construction


[![DOI](https://zenodo.org/badge/362991463.svg)](https://zenodo.org/badge/latestdoi/362991463)

We perform extensive real data analyses to systematically evalute several key factors in supervised celltyping tasks, including: feature selection, prediction method, data preprocessing and the choice of the reference dataset. In our paper, we benchmark **9** classifiers (MLP, SVM with Radial Basis Function kernel , SVM with linear kernel, Random Forest, GEDFN, scmap, CHETAH, ItClust and MARS) along with **6** feature selection strategies. We also investigate **the impact of data preprocessing** (e.g., batch effect correction between reference and target; data imputation on both reference and target), **size of reference**, **number of cell types in reference** and **different annotation methods in target**. Furthermore, we focus on showing how **discrepancies between reference and target data** would affect the prediction performance. In the end, we explore the **strategies of pooling and purifying** reference data. 

Our paper is currently accepted by Genome Biology and will be online soon.


### Repo description

We provide all scripts here for running our prediction and analysis results which can be used to reproduce the results in the paper.

#### Environment setup

In our analysis, we in total set up three conda virtual environments in case some underlying packages are incompatible. 

We have stored the package requirements under each virtual environment to `envs/`and the corresponding environment can be installed by `conda create -n NAME --file xx_env.txt`.  Be aware that R packages installed under the virtual environment will not be listed in the requirement file. We install all R packages (shown below) under the `celltyping` environment:

- [FEAST](https://github.com/suke18/FEAST): Accurate feature selection improves single-cell RNA-seq cell clustering.
-  [scmap](https://bioconductor.org/packages/release/bioc/html/scmap.html): scmap: projection of single-cell RNA-seq data across data sets.
- [CHETAH](http://www.bioconductor.org/packages/release/bioc/html/CHETAH.html): Fast and accurate scRNA-seq cell type identification.
- [Harmony](https://github.com/immunogenomics/harmony): Fast, sensitive and accurate integration of single-cell data with Harmony.
- [fastMNN](http://bioconductor.org/packages/devel/bioc/vignettes/batchelor/inst/doc/correction.html):Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.
- [SAVER](https://github.com/mohuangx/SAVER): SAVER: gene expression recovery for single-cell RNA sequencing.

#### Code structure

- `xx.slurm`: These files are used to run experiments on the cluster. If you are also using a computation cluster, you may change some `#SBATCH` parameters and directly run your own clusters. 
- `pipelines/`: This folder contains all written pipelines along with some auxiliary and analytical scripts. The `xx_pipeline.py` are all written in Python with [argparse](https://docs.python.org/3/library/argparse.html) to parse argument such as classifier name, feature selection method and other parameters. The `xx_analysis.py` scripts are used to analyze the result and plot figures. The `xx_utils.py` are used as auxiliary functions to train the classifiers or load the data. `Rcode` folder contain FEAST, scmap and CHETAH which are also embedded in Python script to execute.
- `preprocessing/`: This folder contain all data preprocessing including loading the data into an [AnnData](https://anndata.readthedocs.io/en/latest/) object and use [scanpy](https://scanpy.readthedocs.io/en/stable/) to preprocess the data.
- `test_xx/`: These folders contain the model implementations.
  - `test_GEDFN`: MLP and [GEDFN](https://github.com/yunchuankong/GEDFN) are implemented here. We have made some minor modifications to GEDFN to adapt to our pipeline.
  - `test_ItClust`: [ItClust](https://github.com/jianhuupenn/ItClust) is included here. We also made some minor modifications to adapt to our pipeline.
  - `test_batcheffect`: Harmony and fastMNN are included here.
  - `test_imputation`: SAVER, [MAGIC](https://github.com/KrishnaswamyLab/MAGIC) and [svVI](https://github.com/YosefLab/scvi-tools) are included here.
  - SVM and Random Forest are directly implemented by [scikit-learn](https://scikit-learn.org/stable/) package.

#### Datasets

- Mouse brain

| Dataset Name | Dataset Description | Protocol | No. Cells | No. Major cell types (subtypes) |
| --- | --- | --- | --- | --- |
| Mouse brain FC <sup>1</sup> | [GSE116470](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116470), Frontal cortex brain region, 7 male adult mice subjects | Drop-seq | 71,639 | 14 (81) |
| Mouse brain HC <sup>1</sup> | [GSE116470](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116470), Hippocampus cortex brain region, 6 male adult mice subjects | Drop-seq | 53,204 | 12 (103) |
| Mouse brain pFC <sup>2</sup> | [GSE124952](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124952), 6 saline-treated adult mice (2 in each 3 timepoints: control, 48h after cocaine withdrawal (CW), 15 days after CW) | 10X Chromium | 11,886 | 8 (9) |  
| Mouse brain cortex <sup>3</sup> | [SCP425](https://singlecell.broadinstitute.org/single_cell/study/SCP425/single-cell-comparison-cortex-data), cortex1 and cortex2 samples from one-month old mice | DroNc-seq | 1,452 (cortex1), 892 (cortex2) | 8 |
| Mouse brian Allen <sup>4</sup> | [NeMO: dat-jb2f34y](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x), 3 male adult mice with frontal cortex extracted | 10X Chromium | 65,944 | 8 (47) |


- Human peripheral blood mononuclear cells (PBMC)

| Dataset Name | Dataset Description | Protocol | No. Cells | No. Major cell types (subtypes) |
| --- | --- | --- | --- | --- |
| Human PBMC lupus <sup>5</sup> | [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583), batch1, 8 SLE patients | 10X Chromium | 12,544 | 6 (8) |
| Human PBMC lupus <sup>5</sup> | [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583), batch2, 8 SLE patients untreated for 6 hours | 10X Chromium | 12,138 | 6 (8) |
| Human PBMC lupus <sup>5</sup> | [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583), batch2, 8 SLE patients activated by IFN-β for 6 hours | 10X Chromium | 12,167 | 6 (8) |
| Human PBMC protocols <sup>3</sup> | [SCP424](https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data), pooled frozen 25million pbmc1 and within 4-hour fresh blood pbmc2 | Smart-seq2/CEL-Seq2/10X Chromium (v2) | 6,814 (pbmc1), 223 (pbmc2) | 6 (9) | 
| Human PBMC FACS <sup>6</sup> | [10X Genomics Datasets](https://www.10xgenomics.com/resources/datasets), fresh healthy Donor A with 10 bead-enriched subpopulations | FACS | 94,655 | 5 (10) |


- Human pancreas

| Dataset Name | Dataset Description | Protocol | No. Cells | No. Major cell types (subtypes) |
| --- | --- | --- | --- | --- |
| Human pancreas Muraro <sup>7</sup> | [GSE85241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241), 4 dead donors (1 female, 3 males; variation in Age and BMI), 8 libraries | CEL-Seq2 | 2,018 | 6 |
| Human pancrease Segerstolpe <sup>8</sup> | [E-MTAB-5061](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/), 6 healthy and 4 T2D individuals (variation in healthy gender and age, BMI) | Smart-Seq2 | 2,038 | 6 |
| Human pancreas Xin <sup>9</sup> | [GSE81608](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81608), 12 Healthy and 6 T2D donors (balanced gender, varied age, BMI, weight) | SMARTer | 1,492 | 6 | 


**References**
1. Saunders A, Macosko EZ, Wysoker A, Goldman M, Krienen FM, Rivera H de, et al. Molecular Diversity and Specializations among the Cells of the Adult Mouse Brain. Cell. 2018;174:1015-1030.e16. 
2. Bhattacherjee A, Djekidel MN, Chen R, Chen W, Tuesta LM, Zhang Y. Cell type-specific transcriptional programs in mouse prefrontal cortex during adolescence and addiction. Nature communications. Nature Publishing Group; 2019;10:1–18. 
3. Ding J, Adiconis X, Simmons SK, Kowalczyk MS, Hession CC, Marjanovic ND, et al. Systematic comparison of single-cell and single-nucleus RNA-sequencing methods. Nature biotechnology. Nature Publishing Group; 2020;38:737–46. 
4. Yao Z, van Velthoven CTJ, Nguyen TN, Goldy J, Sedeno-Cortes AE, Baftizadeh F, et al. A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation. Cell. 2021;184:3222-3241.e26. 
5. Kang HM, Subramaniam M, Targ S, Nguyen M, Maliskova L, McCarthy E, et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nature biotechnology. Nature Publishing Group; 2018;36:89. 
6. Zheng GX, Terry JM, Belgrader P, Ryvkin P, Bent ZW, Wilson R, et al. Massively parallel digital transcriptional profiling of single cells. Nature communications. Nature Publishing Group; 2017;8:1–12. 
7. Muraro MJ, Dharmadhikari G, Grün D, Groen N, Dielen T, Jansen E, et al. A single-cell transcriptome atlas of the human pancreas. Cell systems. Elsevier; 2016;3:385-394. e3. 
8. Segerstolpe Å, Palasantza A, Eliasson P, Andersson E-M, Andréasson A-C, Sun X, et al. Single-cell transcriptome profiling of human pancreatic islets in health and type 2 diabetes. Cell metabolism. Elsevier; 2016;24:593–607. 
9. Xin Y, Kim J, Okamoto H, Ni M, Wei Y, Adler C, et al. RNA sequencing of single human islet cells reveals type 2 diabetes genes. Cell metabolism. Elsevier; 2016;24:608–15. 

Our code is released under MIT license.
