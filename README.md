# GSFA: Bayesian sparse factor analysis on single-cell CRISPR screen data

Single-cell RNA-seq with multiplexed CRISPR screening (e.g. CROP-seq, Perturb-seq) 
is a powerful tool to study the transcriptomic effects of genetic perturbations, 
bringing insights into gene regulation. However, traditional analysis of these data 
poses various statistical and interpretation challenges.

`GSFA` (guided sparse factor analysis) is a methodology (and available as an R package) 
that is particularly suited for single-cell CRISPR screening data. 

Our Bayesian hierarchical model performs factorization of the gene expression matrix with 
a sparse prior on the gene loading to factors. Importantly, the factors in a sample 
depend on the sample-level perturbations through a linear regression model with shrinkage, 
capturing the broad transcriptomic effects of perturbations. In this way, the inference 
of factors are, to some extent, guided by the given perturbations.

Provided with a normalized gene expression matrix and a gRNA design matrix,   
GSAF can identify coordinated factors (gene modules) and estimate their association with 
given genetic perturbations in a joint probabilistic framework. 
Differentially expressed genes under each gRNA perturbation can be detected 
by thresholding the local false sign rate (_lfsr_). The biological function 
of each factor (gene module) can be interpreted through gene ontology/pathway enrichment analysis 
on genes with nonzero weights in the factor.  

![input output](man/figures/schematic.png)

This package is developed by Yifan Zhou from the
[He Lab](http://xinhelab.org) at the University of Chicago.


## Installation

To install the development version of the `GSFA` package from this repository, run:

```R
install.packages("devtools")
devtools::install_github("gradonion/GSFA")
```
