# GSFA: Bayesian sparse factor analysis on single-cell CRISPR screen data

Single-cell RNA-seq with multiplexed CRISPR screening (e.g. CROP-seq, Perturb-seq) 
is a powerful tool to study the transcriptomic effects of genetic perturbations, 
bringing insights into gene regulation. However, traditional analysis of these data 
poses various statistical and interpretation challenges.

`GSFA` (Guided Sparse Factor Analysis) assumes that the perturbation of a target gene 
affects certain latent factors, which in turn changes the expression of individual genes, 
and identifies genetically controlled factors that are associated with the perturbation 
in a joint statistical framework. It also summarizes the effects of a perturbation 
on individual genes as the sum of effects mediated by all the factors.

Provided with a normalized gene expression matrix and a perturbation matrix 
(of which cells contain which type of gRNAs) a single-cell CRISPR screening experiment, 
`GSFA` can   
(1) identify coordinated factors (gene modules) and their associations with the 
genetic perturbations;    
(2) interpret the biological meanings of factors (gene modules) through gene 
ontology/pathway enrichment analysis thanks to the sparse gene weights on factors;    
(3) detect differentially expressed genes under each genetic perturbation by 
thresholding the local false sign rate (_LFSR_).

![input output](man/figures/schematic.png)

## Citing this work

If you find the GSFA package or any of the source code in this
repository useful for your work, please cite:

> Yifan Zhou, Kaixuan Luo, Mengjie Chen and Xin He. 
> A novel Bayesian factor analysis method improves detection of genes and 
> biological processes affected by perturbations in single-cell CRISPR screening. 
> *bioRxiv* doi: 10.1101/2022.02.13.480282 (2021).

## License

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].

## Installation

To install the development version of the `GSFA` package from Github, run:

```R
install.packages("devtools")
devtools::install_github("gradonion/GSFA")
```

Note that installing the package will require a C++ compiler setup that is appropriate for the version of R installed on your computer.

## Using the package

For guidance on using GSFA to analyze single-cell CRISPR screen data, please see the code in this 
[repository][paper_github], and analysis results for the paper [here][result_website].

## Credits

The GSFA package is developed by Yifan Zhou from the
[He Lab](http://xinhelab.org) at the University of Chicago.

[mit-license]: https://opensource.org/licenses/mit-license.html
[paper_github]: https://github.com/gradonion/GSFA_paper/
[result_website]: https://gradonion.github.io/GSFA_paper/
