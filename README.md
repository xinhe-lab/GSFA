# GSFA

GSFA (guided sparse factor analysis) is a Bayesian matrix factorization model 
that can identify coordinated gene modules and estimate their association with 
given genetic perturbations in a joint probabilistic framework. It is particularly
suited for single-cell RNA-seq data with multiplexed CRISPR screening (e.g. CROP-seq, 
Perturb-seq data).

Our hierarchical model performs factorization of the gene expression matrix with 
a sparse prior on the gene loading to factors. Importantly, the factors in a sample 
depend on the sample-level perturbations through a linear regression model with shrinkage, 
capturing the broad transcriptomic effects of perturbations. In this way, the inference 
of factors are, to some extent, guided by the given perturbations.

This package is developed by Yifan Zhou from the
[He Lab](http://xinhelab.org) at the University of Chicago.

## Quick Start

To automatically retrieve and install `GSFA` from this repository, run:

```R
devtools::install_github("gradonion/GSFA")
```
