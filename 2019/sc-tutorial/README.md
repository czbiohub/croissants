# Scripts for "Current best-practices in single-cell RNA-seq: a tutorial"

![image](https://drive.google.com/uc?export=view&id=1YoX3F8gNGH5K0AFu4wFdd5ccEXT-mP_J)

This repository is complementary to the publication:

M.D. Luecken, F.J. Theis, "Current best practices in single-cell RNA-seq analysis: a tutorial", Molecular Systems Biology 15(6) (2019): e8746

The paper was recommended on F1000 prime as being of special significance in the field.

<a href="https://f1000.com/prime/736010853?bd=1 &ui=150648 " target="_blank"><img src="https://s3.us-east-1.amazonaws.com/cdn.f1000.com/images/badges/badgef1000.gif" alt="Access the recommendation on F1000Prime" id="bg" /></a>


The repository contains:
- scripts to generate the paper figures
- a case study which complements the manuscript
- the code for the marker gene detection study from the supplementary material

The main part of this repository is a case study where the best-practices established in the manuscript are applied to a mouse intestinal epithelium regions dataset from Haber et al., Nature 551 (2018) available from the GEO under GSE92332. This case study can be found in different versions in the `latest_notebook/` and `old_releases/` directories.

The scripts in the plotting_scripts/ folder reproduce the figures that are shown in the manuscript and the supplementary materials. These scripts contain comments to explain each step. Each figure that does not have a corresponding script in the plotting_scripts/ folder was taken from the case study or the marker gene study.

In case of questions or issues, please get in touch by posting an issue in this repository.

If the materials in this repo are of use to you, please consider citing the above publication.


## Environment set up

To run the tutorial case study, several packages must be installed. As both `R` and `python` packages are required. To facilitate setup, there is a `requirements.txt` file, which contains all pip installable dependencies. R dependencies, which are not already available, must be installed in the notebook.


To install the `python` packages, run the following command.

```
pip3 install -r requirements.txt
```

To install the `R` dependencies run the install commands already in the notebook:

```
install.packages(c('devtools', 'gam', 'RColorBrewer', 'BiocManager'))
update.packages(ask=F)
BiocManager::install(c("scran","MAST","monocle","ComplexHeatmap","slingshot"), version = "3.9")
```
 
## Downloading the data

As mentioned above the data for the case study comes from GSE92332. To run the case study as shown, you must download this data and place it in the correct folder. Unpacking the data requires `tar` and `gunzip`, which should already be available on most systems. If you are cloning the github repository and have the case study script in a `latest_notebook/` folder, then from the location where you store the case study ipynb file, this can be done via the following commands:

```
cd ../  #To get to the main github repo folder
mkdir -p data/Haber-et-al_mouse-intestinal-epithelium/
cd data/Haber-et-al_mouse-intestinal-epithelium/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_RAW.tar
mkdir GSE92332_RAW
tar -C GSE92332_RAW -xvf GSE92332_RAW.tar
gunzip GSE92332_RAW/*_Regional_*
```

The annotated dataset with which we briefly compare the results at the end of the notebook, is available from the same GEO accession ID (GSE92332). It can be obtained using the following command:

```
cd data/Haber-et-al_mouse-intestinal-epithelium/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_Regional_UMIcounts.txt.gz
gunzip GSE92332_Regional_UMIcounts.txt.gz
```


## Case study notes

We have noticed that results such as visualization, dimensionality reduction, and clustering (and hence all downstream results as well) can give slightly different results on different systems. This has to do with the numerical libraries that are used in the backend. Thus, we cannot guarantee that a rerun of the notebook will generate exactly the same clusters.

While all results are qualitatively similar, the assignment of cells to clusters especialy for stem cells, TA cells, and enterocyte progenitors can differ between runs across systems. To show the diversity that can be expected, we have uploaded shortened case study notebooks to the `alternative_clustering_results/` folder.

Note that running `sc.pp.pca()` with the parameter `svd_solver='arpack'` drastically reduces the variability between systems, however the output is not exactly the same.


## Adapting the pipeline for other datasets:

The pipeline was designed to be easily adaptable to new datasets. However, there are several limitations to the general applicability of the current workflow. When adapting the pipeline for your own dataset please take into account the following:

1. Sparse data formats are not supported by `rpy2` and therefore do not work with any of the integrated R commands. Datasets can be turned into a dense format using the code: `adata.X = adata.X.toarray()`

2. The case study assumes that the input data is count data obtained from a single-cell protocol with UMIs. If the input data is full-length read data, then one could consider replacing the normalization method with another method that includes gene length normalization (e.g., TPM).


## Manual installation of package requirements

The following packages are required to run the first version of the case study notebook. For further versions see the README.md in the latest_notebook/ and old_releases/ folders.

General:
- Jupyter notebook
- IRKernel
- rpy2
- R >= 3.4.3
- Python >= 3.5

Python:
- scanpy==1.4.3
- anndata==0.6.21
- numpy
- scipy
- pandas
- seaborn
- louvain>=0.6
- python-igraph==0.7.1.post4
- gprofiler-official (from Case study notebook 1906 version)
- python-gprofiler from Valentine Svensson's github (vals/python-gprofiler)
  - only needed for notebooks before version 1906
- ComBat python implementation from Maren Buettner's github (mbuttner/maren_codes/combat.py)
  - only needed for scanpy versions before 1.3.8 which don't include `sc.pp.combat()`

R:
- scater
- scran
- MAST
- gam
- slingshot (change DESCRIPTION file for R version 3.4.3)
- monocle 2
- limma
- ComplexHeatmap
- RColorBrewer
- clusterExperiment
- ggplot2

## Possible sources of error in the manual installation:

#### For R 3.4.3:

When using Slingshot in R 3.4.3, you must pull a local copy of slingshot via the github repository and change the `DESCRIPTION` file to say `R>=3.4.3` instead of `R>=3.5.0`.

#### For R >= 3.5 and bioconductor >= 3.7:

The clusterExperiment version that comes for bioconductor 3.7 has slightly changed naming convention. `clusterExperiment()` is now called `ClusterExperiment()`. The latest version of the notebook includes this change, but when using the original notebook, please note that this may throw an error.

#### For rpy2 < 3.0.0:

Pandas 0.24.0 is not compatible with rpy2 < 3.0.0. When using old versions of rpy2, please downgrade pandas to 0.23.X. Please also note that Pandas 0.24.0 requires anndata version 0.6.18 and scanpy version > 1.37.0.

#### For enrichment analysis with g:profiler:

Ensure that the correct g:profiler package is used for the notebook. Notebooks until 1904 use `python-gprofiler` from valentine svensson's github, and Notebooks from 1906 use the `gprofiler-official` package from the g:profiler team.
