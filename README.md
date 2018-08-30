# `corncob` QIIME2 plugin

[![Build Status](https://travis-ci.org/bryandmartin/q2-corncob.svg?branch=master)](https://travis-ci.org/bryandmartin/q2-corncob)

This repository contains the `corncob` QIIME 2 plugin. `corncob` is in active development and is available in `R` (https://github.com/bryandmartin/corncob) or as a QIIME 2 plugin (q2-corncob).


`corncob` is based in `R` and requires installation of dependencies `VGAM`, `devtools`, `magrittr`, `phyloseq`, and `dplyr` into your `conda` environment before installing `corncob`. Please refer to the following instructions on how to install `corncob` and its dependencies.


### Activate your QIIME Environment

- Here we activate our example version of QIIME, `qiime2-2018.6`. If you're not sure what your current version of QIIME is you can run `conda env list` in the command line to see a list of installed QIIME environments.

```
source activate qiime2-2018.6
```


### Install `corncob` dependencies

(Expected installation time ~3-5 minutes)

```
conda install -c bioconda -c conda-forge bioconductor-phyloseq r-devtools r-magrittr r-dplyr r-vgam unzip
```

- Note: When installing select `y` to proceed with installation when prompted.


### Install  `corncob`  and `q2-corncob`

```
pip install git+https://github.com/statdivlab/q2-corncob.git
qiime dev refresh-cache
```


### Check that `corncob` is installed <br>

```
qiime corncob --help
```


# QIIME 2 Tutorial: Using `q2-corncob`
This is a Community Tutorial for q2-corncob within the qiime2-2018.6 release.

`corncob` is an individual taxon regression model that uses abundance tables and sample data. `corncob` is able to model differential abundance and differential variability and address the following statistical challenges with modeling micriobial relative abundance:
- different sequencing depth
- excessive zeros from unobserved taxa
- high variability of empirical relative abundance (overdispersion)
- within-taxon correlation
- hypothesis testing with categorical and continuous covariates

A vignette on how to use `corncob` in `R` can be found [here]().

Currently, q2-corncob has made available the following functions within `corncob`:
- `differentialtest()`

### Citing `corncob`
Please cite the following when using `corncob`...
- `differentialtest()`

### How to use `q2-corncob`
For this tutorial we will be using data from the "Moving Pictures" tutorial. q2-corncob requires input of a FeatureTable, Metadata, Taxonomy, and a covariate of interest.

<a href="https://github.com/statdivlab/q2-corncob/blob/master/data/table.qza?raw=true" download>table.qza</a>
<br><a href="https://github.com/statdivlab/q2-corncob/blob/master/data/taxonomy.qza?raw=true" download>taxonomy.qza</a>
<br><a href="https://raw.githubusercontent.com/statdivlab/q2-corncob/master/data/metadata.tsv" download>metadata.tsv</a>


Let's say that we are interested in seeing if there are ASV's that are differentially abundant or differentially variable across groups of ReportedAntibioticUsage.
```
qiime corncob differentialtest \
--i-table table.qza \
--m-metadata-file metadata.tsv \
--p-variable ReportedAntibioticUsage \
--i-taxonomy taxonomy.qza \
--o-output corncobresults
```
Our results show a table of features, taxonomic assignment, and fdr controlled p-values for differential abundance and differential variance.

<a href="https://github.com/statdivlab/q2-corncob/blob/master/data/taxonomy.tsv" download>results</a>
