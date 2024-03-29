% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EPIC.R
\name{process_scRNA}
\alias{process_scRNA}
\title{Transcriptomic data processing for scRNA-seq}
\usage{
process_scRNA(SeuratObject, meta_ct_col)
}
\arguments{
\item{SeuratObject}{a Seurat object from a gene expression matrix, returned from \code{CreateSeuratObject}.}

\item{meta_ct_col}{column name of cell types in the meta data matrix.}
}
\value{
A list with components
 \item{scRNA.rpkm}{A data frame of cell-type-specific RPKMs with log2 transformation}
 \item{scRNA.loc}{A data frame of feature dictionary for top 8000 highly variable genes}
}
\description{
For scRNA-seq data, we follow the Seurat pipeline to perform gene- and cell-wise
quality controls and focus on the top 8000 highly variable genes.
Cell-type-specific RPKMs are calculated by combining reads or UMI counts from
all cells pertaining to a specific cell type, followed by log2 transformation
with an added pseudo-count.
}
\examples{
\dontrun{
 scRNA.object <- process_scRNA(SeuratObject, meta_ct_col = "cell.type")
 scRNA.rpkm <- scRNA.object$scRNA.rpkm
 scRNA.loc <- scRNA.object$scRNA.loc
}

}
\author{
Rujin Wang \email{rujin@email.unc.edu}
}
