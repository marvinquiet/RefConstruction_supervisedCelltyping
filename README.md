# Reference construction strategies for single-cell supervised celltyping
We perform extensive real data analyses to systematically evalute several key factors in supervised celltyping tasks, including: feature selection, prediction method, data preprocessing and the choice of the reference dataset. In our paper, we benchmark **9** classifiers (MLP, SVM with Radial Basis Function kernel , SVM with linear kernel, Random Forest, GEDFN, scmap, CHETAH, ItClust and MARS) along with **6** feature selection strategies. We also investigate **the impact of data preprocessing** (e.g., batch effect correction between reference and target; data imputation on both reference and target), **size of reference**, **number of cell types in reference** and **different annotation methods in target**. Furthermore, we focus on showing how **discrepancies between reference and target data** would affect the prediction performance. In the end, we explore the **strategies of pooling and purifying** reference data. 

Our paper is currently under review.



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

Once our paper gets published, we will update this section