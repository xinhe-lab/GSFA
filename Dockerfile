FROM rocker/tidyverse:4.0.4
MAINTAINER Yifan Zhou, zhouyf@uchicago.edu

RUN apt-get update && apt-get install -y --no-install-recommends \
        libgit2-dev \
        libxml2-dev \
        curl \
        tree \
        jq \
        htop \
        texinfo \
        vim \
        man-db \
        less \
        g++-11

ENV MRAN_BUILD_DATE=2020-09-01
# Install Basic Utility R Packages
RUN install2.r -r https://cran.microsoft.com/snapshot/${MRAN_BUILD_DATE} \
    --error \
    Rcpp \
    RcppArmadillo \
    data.table \
    magrittr \
    reshape2 \
    ggplot2 \
    knitr \
    rmarkdown \
    devtools

RUN installGithub.r gradonion/GSFA
