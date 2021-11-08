rm(list = ls())
curr_dir = "/pine/scr/r/u/rujin/EPIC"
input_dir = "/pine/scr/r/u/rujin/scRNAseqGWAS"
func_dir = "/pine/scr/r/u/rujin/myfunction/scRNAseqGWAS"


setwd(input_dir)
library(data.table)
library(dplyr)
library(GWASTools)
library(Matrix)
library(MASS)
library(psych)
library(ggplot2)
library(ggpubr)
library(scales)

source(file = file.path(func_dir, "epic.R"))

type = "POET_QC"

trait = "blah_trait_blah"
maf_thres = 0.01
sliding_step = 1
sliding_size = 10
Cmin = 1.5

out_dir <- file.path(curr_dir, "gtex_v8", "my_x_results", trait, paste0(type, "_", sliding_size, "_", Cmin, "_results"))
in_dir <- file.path(curr_dir, "gtex_v8", "my_x_results", trait, paste0(type, "_", sliding_size, "_", Cmin))
if(!dir.exists(out_dir)){
  dir.create(out_dir, recursive = TRUE)
}
if(!dir.exists(in_dir)){
  dir.create(in_dir, recursive = TRUE)
}

loc <- read.table(file = file.path(curr_dir, "gtex_v8", paste0("gtex_v8.gene.noMHC.loc")), stringsAsFactors = FALSE)
head(loc)
dim(loc)
head(loc, 20)

# sort the genes by genomic positions!!!
loc <- loc[order(loc$V2, loc$V3, loc$V4), ]
head(loc)
dim(loc)
tail(loc, 20)


load(file = file.path(curr_dir, "trait_results", trait, paste0(trait, "_common_", maf_thres, "_chisq_spd_QC.rda")))
prunein_gene_pval <- prunein_gene_pval[which(!is.na(match(names(prunein_gene_pval), loc$V1)))]
length(prunein_gene_pval)

loc <- loc[which(!is.na(match(loc$V1, names(prunein_gene_pval)))), ]
dim(loc)

# reorder gene_pval list by position!!!!!!!!!!!!
prunein_gene_pval <- prunein_gene_pval[match(loc$V1, names(prunein_gene_pval))]
all(names(prunein_gene_pval)==loc$V1)

prunein_POET_gene_pval <- prunein_gene_pval
prunein.pval <- sapply(prunein_gene_pval, function(z){z$pval}, USE.NAMES = FALSE)
rm("prunein_gene_pval")

prunein_genes = names(prunein_POET_gene_pval)
length(prunein_genes)
head(prunein_genes)

prunein_gene_corr = vector("list", length = length(prunein_POET_gene_pval))
names(prunein_gene_corr) = names(prunein_POET_gene_pval)

for(chr in 1:22){
  nsnps <- sum(loc$V2==chr)
  genes.list <- loc$V1[which(loc$V2==chr)]
  loc.temp <- loc[which(loc$V2==chr), ]
  for (i in seq(1, (nsnps - sliding_size + 1), sliding_step)) {
    load(file = file.path(in_dir, paste0(trait, "_", maf_thres, "_LD_chr", chr, "_", genes.list[i], "_", type, ".rda")))
    sw.gene <- c(read.table(file.path(in_dir, paste0(trait, "_", maf_thres, "_LD_chr", chr, "_", genes.list[i], "_", type, ".txt")), 
                            stringsAsFactors = FALSE)$V1)
    gene.name = genes.list[i]
    
    cat(i, ": ", gene.name, " chr = ", chr, "\n")
    
    snps = prunein_POET_gene_pval[[gene.name]]$snps
    prunein_POET_gene_pval[[gene.name]][["LD"]] <- corr.POET[snps, snps, drop = FALSE]
    LD_inverse = compute_ld_inverse(prunein_POET_gene_pval[[gene.name]][["LD"]], tol = 1e-5)
    prunein_POET_gene_pval[[gene.name]][["inv_LD"]] = LD_inverse$inv_mat
    prunein_POET_gene_pval[[gene.name]][["eigenvalue"]] = LD_inverse$ev
    prunein_POET_gene_pval[[gene.name]][["df"]] = LD_inverse$df
    z_mat = as.matrix(prunein_POET_gene_pval[[gene.name]][["zhat"]])
    prunein_POET_gene_pval[[gene.name]][["chisq"]] = as.numeric(t(z_mat) %*% prunein_POET_gene_pval[[gene.name]][["inv_LD"]] %*% z_mat)
    prunein_POET_gene_pval[[gene.name]][["pval"]] = exp(pchisq(prunein_POET_gene_pval[[gene.name]][["chisq"]], df = prunein_POET_gene_pval[[gene.name]][["df"]], log.p = TRUE, lower.tail = FALSE))
    
    # gene-gene correlation
    gene1 = gene.name
    gene1.snps <- prunein_POET_gene_pval[[gene1]][["snps"]]
    
    st <- loc$V3[as.character(loc$V1)==gene1]
    ed <- loc$V4[as.character(loc$V1)==gene1]
    
    temp <- loc.temp[-c(1:i), , drop = FALSE]
    close.gene <- as.character(temp$V1[temp$V2 == chr & temp$V3 < ed + 5000000 & temp$V4 > st - 5000000]) # 5MB
    close.gene <- close.gene[close.gene %in% sw.gene]
    close.gene <- close.gene[close.gene %in% names(prunein_POET_gene_pval)]
    
    if(length(close.gene)>0){
      
      prunein_gene_corr[[gene1]] = vector("list", length = length(close.gene))
      names(prunein_gene_corr[[gene1]]) = close.gene
      
      if(length(close.gene)>0){
        for(gene2 in close.gene){
          gene2.snps <- prunein_POET_gene_pval[[gene2]][["snps"]]
          R <- corr.POET[c(gene1.snps, gene2.snps), c(gene1.snps, gene2.snps), drop = FALSE]
          
          p <- length(prunein_POET_gene_pval[[gene1]]$snps)
          q <- length(prunein_POET_gene_pval[[gene2]]$snps)
          gene1.invR.sqrt = compute_ld_inverse_sqrt(R[1:p, 1:p, drop = FALSE], tol = 1e-5)
          gene2.invR.sqrt = compute_ld_inverse_sqrt(R[(p+1):(p+q), (p+1):(p+q), drop = FALSE], tol = 1e-5)
          
          three.parts = gene1.invR.sqrt %*% R[1:p, (p+1):(p+q)] %*% gene2.invR.sqrt
          R.star <- diag(p+q)
          R.star[1:p, (p+1):ncol(R.star)] <- three.parts
          R.star[(p+1):ncol(R.star), 1:p] <- t(three.parts)
          
          R.star.chol <- cholesky_decomposition(R.star)
          cov.chisq <- compute_covariance_cholesky(Lmat = R.star.chol, p = p, q = q)
          corr.chisq <- cov.chisq / (sqrt(2*p) * sqrt(2*q))
          
          prunein_gene_corr[[gene1]][[gene2]]$R <- R
          prunein_gene_corr[[gene1]][[gene2]]$cov.chisq <- cov.chisq
          prunein_gene_corr[[gene1]][[gene2]]$corr.chisq <- corr.chisq
          
          cat("\t corr with ", gene2, "=", corr.chisq, "\n")
        }
      }
    }
  }
  
  last20 <- length(genes.list) - sliding_size + 2
  load(file = file.path(in_dir, paste0(trait, "_", maf_thres, "_LD_chr", chr, "_", genes.list[(last20 - 1)], "_", type, ".rda")))
  for (i in seq(last20, length(genes.list), sliding_step)) {
    gene.name = genes.list[i]
    cat(i, ": ", gene.name, " chr = ", chr, "\n")
    
    snps = prunein_POET_gene_pval[[gene.name]]$snps
    prunein_POET_gene_pval[[gene.name]][["LD"]] <- corr.POET[snps, snps, drop = FALSE]
    LD_inverse = compute_ld_inverse(prunein_POET_gene_pval[[gene.name]][["LD"]], tol = 1e-5)
    prunein_POET_gene_pval[[gene.name]][["inv_LD"]] = LD_inverse$inv_mat
    prunein_POET_gene_pval[[gene.name]][["eigenvalue"]] = LD_inverse$ev
    prunein_POET_gene_pval[[gene.name]][["df"]] = LD_inverse$df
    z_mat = as.matrix(prunein_POET_gene_pval[[gene.name]][["zhat"]])
    prunein_POET_gene_pval[[gene.name]][["chisq"]] = as.numeric(t(z_mat) %*% prunein_POET_gene_pval[[gene.name]][["inv_LD"]] %*% z_mat)
    prunein_POET_gene_pval[[gene.name]][["pval"]] = exp(pchisq(prunein_POET_gene_pval[[gene.name]][["chisq"]], df = prunein_POET_gene_pval[[gene.name]][["df"]], log.p = TRUE, lower.tail = FALSE))
    
    # gene-gene correlation
    gene1 = gene.name
    gene1.snps <- prunein_POET_gene_pval[[gene1]][["snps"]]
    
    st <- loc$V3[as.character(loc$V1)==gene1]
    ed <- loc$V4[as.character(loc$V1)==gene1]
    
    temp <- loc.temp[-c(1:i), , drop = FALSE]
    close.gene <- as.character(temp$V1[temp$V2 == chr & temp$V3 < ed + 5000000 & temp$V4 > st - 5000000]) # 5MB
    close.gene <- close.gene[close.gene %in% sw.gene]
    close.gene <- close.gene[close.gene %in% names(prunein_POET_gene_pval)]
    
    if(length(close.gene)>0){
      
      prunein_gene_corr[[gene1]] = vector("list", length = length(close.gene))
      names(prunein_gene_corr[[gene1]]) = close.gene
      
      if(length(close.gene)>0){
        for(gene2 in close.gene){
          gene2.snps <- prunein_POET_gene_pval[[gene2]][["snps"]]
          R <- corr.POET[c(gene1.snps, gene2.snps), c(gene1.snps, gene2.snps), drop = FALSE]
          
          p <- length(prunein_POET_gene_pval[[gene1]]$snps)
          q <- length(prunein_POET_gene_pval[[gene2]]$snps)
          gene1.invR.sqrt = compute_ld_inverse_sqrt(R[1:p, 1:p, drop = FALSE], tol = 1e-5)
          gene2.invR.sqrt = compute_ld_inverse_sqrt(R[(p+1):(p+q), (p+1):(p+q), drop = FALSE], tol = 1e-5)
          
          three.parts = gene1.invR.sqrt %*% R[1:p, (p+1):(p+q)] %*% gene2.invR.sqrt
          R.star <- diag(p+q)
          R.star[1:p, (p+1):ncol(R.star)] <- three.parts
          R.star[(p+1):ncol(R.star), 1:p] <- t(three.parts)
          
          R.star.chol <- cholesky_decomposition(R.star)
          cov.chisq <- compute_covariance_cholesky(Lmat = R.star.chol, p = p, q = q)
          corr.chisq <- cov.chisq / (sqrt(2*p) * sqrt(2*q))
          
          prunein_gene_corr[[gene1]][[gene2]]$R <- R
          prunein_gene_corr[[gene1]][[gene2]]$cov.chisq <- cov.chisq
          prunein_gene_corr[[gene1]][[gene2]]$corr.chisq <- corr.chisq
          
          cat("\t corr with ", gene2, "=", corr.chisq, "\n")
        }
      }
    }
  }
}


prunein.pval_shrink <- sapply(prunein_POET_gene_pval, function(z){z$pval}, USE.NAMES = FALSE)
length(prunein.pval)
length(prunein.pval_shrink)


save(prunein_POET_gene_pval, file = file.path(out_dir, paste0(trait, "_common_", maf_thres, "_chisq_", type, ".noMHC.rda")))
save(prunein_gene_corr, file = file.path(out_dir, paste0(trait, "_common_", maf_thres, "_", type, "_gene_corr_5MB.noMHC.rda")))

length(prunein_gene_corr)
prunein_genes = names(prunein_gene_corr)

# read in gene-gene correlation based on chisq statistics
prunein_gene_corr.mat = diag(length(prunein_gene_corr))
rownames(prunein_gene_corr.mat) <- prunein_genes
colnames(prunein_gene_corr.mat) <- prunein_genes
for (i in 1:length(prunein_gene_corr)) {
  if(i%%1000==1){
    cat(i, "\t")
  }
  gene1 = names(prunein_gene_corr)[i]
  if(!is.null(prunein_gene_corr[[gene1]])){
    for (gene2 in names(prunein_gene_corr[[gene1]])) {
      prunein_gene_corr.mat[gene1, gene2] = prunein_gene_corr[[gene1]][[gene2]]$corr.chisq
      prunein_gene_corr.mat[gene2, gene1] = prunein_gene_corr[[gene1]][[gene2]]$corr.chisq
    }
  }
}

save(prunein_gene_corr.mat, file = file.path(out_dir, paste0(trait, "_common_", maf_thres, "_", type, "_gene_corr_mat_5MB.noMHC.rda")))





q("no")


