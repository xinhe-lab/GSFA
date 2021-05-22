# GSFA: Bayesian sparse factor analysis on single-cell CRISPR screen data

Single-cell RNA-seq with multiplexed CRISPR screening (e.g. CROP-seq, Perturb-seq) 
is a powerful tool to study the transcriptomic effects of genetic perturbations, 
bringing insights into gene regulation. However, traditional analysis of these data 
poses various statistical and interpretation challenges.

`GSFA` (guided sparse factor analysis) is a methodology (and available in the form of R package) 
that is particularly suited for single-cell CRISPR screening data. 
It is a Bayesian matrix factorization model 
that can identify coordinated gene modules and estimate their association with 
given genetic perturbations in a joint probabilistic framework.

Our hierarchical model performs factorization of the gene expression matrix with 
a sparse prior on the gene loading to factors. Importantly, the factors in a sample 
depend on the sample-level perturbations through a linear regression model with shrinkage, 
capturing the broad transcriptomic effects of perturbations. In this way, the inference 
of factors are, to some extent, guided by the given perturbations.

<img src="man/figures/schematic.png" align="middle" alt="" width="400" />

This package is developed by Yifan Zhou from the
[He Lab](http://xinhelab.org) at the University of Chicago.


## Installation

To install the development version of the `GSFA` package from this repository, run:

```R
install.packages("devtools")
devtools::install_github("gradonion/GSFA")
```
