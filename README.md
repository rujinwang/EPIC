# EPIC
Inferring relevant tissues and cell types for complex traits in genome-wide association studies


## Authors
Rujin Wang, Danyu Lin, and Yuchao Jiang


## Maintainer
Rujin Wang <rujin@email.unc.edu>


## Installation
From GitHub
```
install.packages('devtools')
devtools::install_github("rujinwang/EPIC/package")
```

## Description
More than a decade of genome-wide association studies (GWAS) have yielded genetic risk variants
that are significantly associated with complex traits. Emerging evidence suggests that the function of 
trait-associated variants likely acts in a tissue- or cell-type- specific fashion. Yet, it remains challenging
to prioritize trait-relevant tissues or cell types to elucidate disease etiology. We present EPIC (cEll tyPe enrIChment), 
a statistical framework that relates large-scale GWAS summary statistics to cell-type-specific omics measurements from
single-cell sequencing. We derive powerful gene-level test statistics for common and rare variants, separately and jointly,
and adopt generalized least squares to prioritize trait-relevant tissues or cell types, while accounting for the correlation
structures both within and between genes.

## Manuscript
Wang R, Lin D, Jiang Y. EPIC: inferring relevant tissues and cell types for complex traits by integrating genome-wide association studies and single-cell RNA sequencing. ***bioRxiv***, 2021. ([link](https://www.biorxiv.org/content/10.1101/2021.06.09.447805v1))

## Vignettes
[HTML](http://htmlpreview.github.io/?https://github.com/rujinwang/EPIC/blob/master/vignettes/EPIC_vignette.html)

