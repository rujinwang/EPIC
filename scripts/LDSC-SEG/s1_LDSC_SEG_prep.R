rm(list = ls())
curr_dir = "/pine/scr/r/u/rujin/EPIC"
input_dir = "/pine/scr/r/u/rujin/scRNAseqGWAS"
func_dir = "/pine/scr/r/u/rujin/myfunction/scRNAseqGWAS"

setwd(input_dir)
library(data.table)
library(rolypoly)
library(dplyr)
library(ggplot2)
library(Matrix)
library(parallel)

source(file = file.path(func_dir, "epic.R"))

ldsc_out_dir <- file.path(curr_dir, "gtex_v8", "ldsc_results")
if(!dir.exists(ldsc_out_dir)){
  dir.create(ldsc_out_dir, recursive = TRUE)
}
if(!dir.exists(file.path(ldsc_out_dir, "SEG"))){
  dir.create(file.path(ldsc_out_dir, "SEG"), recursive = TRUE)
}


gene.expr = read.table(file = file.path(curr_dir, "gtex_v8", paste0("log2_tpm_median.txt")), header = TRUE, row.names = 1)
head(gene.expr); dim(gene.expr)

loc <- read.table(file = file.path(curr_dir, "gtex_v8", paste0("gtex_v8.gene.noMHC.loc")), stringsAsFactors = FALSE)
head(loc)
dim(loc)
head(loc, 20)

gene.expr <- gene.expr[-which(!(rownames(gene.expr) %in% loc$V1)),]
head(gene.expr); dim(gene.expr)

# remove the average columns
gene.expr <- gene.expr[, !(colnames(gene.expr) %in% "Average")]
head(gene.expr); dim(gene.expr)

tstat.mat = compute_tstat(tissues = colnames(gene.expr), gene_mat = gene.expr)
dim(tstat.mat); head(tstat.mat)

# keep top10% SEG
n_genes <- nrow(gene.expr)
n_genes_to_keep <- (n_genes * 0.1) %>% round()

top10pct = apply(t(tstat.mat), 1, function(x) order(-x)[1:n_genes_to_keep])
dim(top10pct)
head(top10pct)
allgene = rownames(gene.expr)
length(allgene)
head(allgene)

for(ct in colnames(gene.expr)){
  gene_top10pct = allgene[top10pct[,ct]]
  write.table(gene_top10pct, file = file.path(ldsc_out_dir, "SEG", paste0("top10pct_GTEx_v8_",ct,".txt")), col.names=F, row.names = F, quote = F)
}


control.genes <- allgene
write.table(control.genes, file = file.path(ldsc_out_dir, "SEG", paste0("control.txt")), col.names=F, row.names = F, quote = F)





q("no")











