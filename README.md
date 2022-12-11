# EPIC
Inferring relevant tissues and cell types for complex traits in genome-wide association studies


## Authors
Rujin Wang, Danyu Lin, and Yuchao Jiang


## Maintainer
Rujin Wang <rujinwang1124@gmail.com>


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
Wang R, Lin D, Jiang Y. EPIC: inferring relevant tissues and cell types for complex traits by integrating genome-wide association studies and single-cell RNA sequencing. ***PLOS Genetics***, 2022. ([link](https://journals.plos.org/plosgenetics/article/authors?id=10.1371/journal.pgen.1010251))

## Vignettes
[HTML](http://htmlpreview.github.io/?https://github.com/rujinwang/EPIC/blob/master/vignettes/EPIC_vignette.html)

