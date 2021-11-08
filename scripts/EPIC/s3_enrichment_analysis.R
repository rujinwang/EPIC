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

load(file = file.path(out_dir, paste0(trait, "_common_", maf_thres, "_chisq_", type, ".noMHC.rda")))
load(file = file.path(out_dir, paste0(trait, "_common_", maf_thres, "_", type, "_gene_corr_mat_5MB.noMHC.rda")))

loc <- read.table(file = file.path(curr_dir, "gtex_v8", paste0("gtex_v8.gene.noMHC.loc")), stringsAsFactors = FALSE)
head(loc)
dim(loc)
head(loc, 20)

# sort the genes by genomic positions!!!
loc <- loc[order(loc$V2, loc$V3, loc$V4), ]
head(loc)
dim(loc)
tail(loc, 20)

# union of genes from rare + common variants
prunein.genes = names(prunein_POET_gene_pval)
length(prunein.genes)

# read in gene expression log input
gene.expr = read.table(file = file.path(curr_dir, "gtex_v8", paste0("log2_tpm_median.txt")), header = TRUE, row.names = 1)
head(gene.expr); dim(gene.expr)
gene.expr.keep <- gene.expr[match(prunein.genes, rownames(gene.expr)),]
head(gene.expr.keep); dim(gene.expr.keep)

all(rownames(gene.expr.keep)==names(prunein_POET_gene_pval))

chisq <- sapply(prunein_POET_gene_pval, function(z){z$chisq}, USE.NAMES = FALSE)
df <- sapply(prunein_POET_gene_pval, function(z){z$df}, USE.NAMES = FALSE)
prunein.pval <- sapply(prunein_POET_gene_pval, function(z){z$pval}, USE.NAMES = FALSE)
length(chisq); length(df); length(prunein.pval)
prunein.pval.zero = which(prunein.pval<1e-300)
chisq[prunein.pval.zero] <- qchisq(1e-300, df = df[prunein.pval.zero], lower.tail = FALSE)
y <- chisq/df
head(y); length(y); summary(y)

# remove the "Average" column
tissues = colnames(gene.expr.keep)[-ncol(gene.expr.keep)]

# Read in gene-gene correlation based on chisq statistics
all(names(prunein_POET_gene_pval)==colnames(prunein_gene_corr.mat))

# Construct a block-diagonal matrix
corr.list <- vector('list', 22)
chol.inv.list <- vector('list', 22)
y.list <- vector('list', 22)
y.reg.list <- vector('list', 22)

y.pval.list <- vector('list', 22)
y.pval.reg.list <- vector('list', 22)

for (chr in 1:22) {
  chr.genes <- loc$V1[loc$V2 == chr]
  block.idx <- match(chr.genes, prunein.genes)[which(!is.na(match(chr.genes, prunein.genes)))]
  corr.list[[chr]] <- prunein_gene_corr.mat[block.idx, block.idx]
  y.list[[chr]] <- y[block.idx]
  chol.inv.list[[chr]] = compute_ld_inverse_sqrt(corr.list[[chr]], tol = 1e-3)
  y.reg.list[[chr]] <- chol.inv.list[[chr]] %*% y.list[[chr]]
  cat("chr", chr, "with corr =", cor(y.list[[chr]], y.reg.list[[chr]]), '\n')
}


y.reg <- unlist(y.reg.list[1:22])
gene.expr.keep.reg <- bdiag(chol.inv.list[1:22]) %*% as.matrix(gene.expr.keep)
gene.expr.keep.reg <- as.matrix(gene.expr.keep.reg)


gene_ensemble <- read.delim(file = file.path(curr_dir, "ENSG.genes.txt"), stringsAsFactors = FALSE, header = TRUE, sep = "\t")
dim(gene_ensemble)
head(gene_ensemble)

chisq.pvals0 <- rep(NA, length(tissues))
names(chisq.pvals0) <- tissues
chisq.pvals <- rep(NA, length(tissues))
names(chisq.pvals) <- tissues

# influential points
pdf(file = file.path(out_dir, paste0(trait, "_", maf_thres, "_common_lm_diagnostics.pdf")), width = 12, height = 12)
for (i in tissues) {
  par(mfrow=c(2,2))
  j = 1
  
  genes <- names(prunein_POET_gene_pval)
  y.qc <- y.reg
  x1.qc <- gene.expr.keep.reg[,i]
  x2.qc <- gene.expr.keep.reg[,"Average"]
  df.qc <- df
  
  cat(i, "\n")
  
  while(j <= 5){
    lm0 = lm(y.qc ~ x1.qc + x2.qc, weights = df.qc/2)
    chisq.pvals0[i] = exp(pt(summary(lm0)$coefficients[2,3], df = summary(lm0)$df[2], lower.tail = FALSE, log.p = TRUE))
    cooksd <- cooks.distance(lm0)
    names(cooksd) <- gene_ensemble$external_gene_name[match(genes, gene_ensemble$ensembl_gene_id)]
    plot(cooksd, main = "Influential Obs by Cook's distance", ylab = paste0("Cook's Distance for ", i))
    cooksd_thres <- 0.2
    abline(h = cooksd_thres, lty = 2, col = "red")
    text(x = 1:length(cooksd), y = cooksd, labels = ifelse(cooksd > cooksd_thres, names(cooksd), ""), col = "red")
    
    # # zoom in cook's distance
    # plot(cooksd, main = "Influential Obs by Cook's distance", ylab = paste0("Cook's Distance for ", i), ylim = c(0, cooksd_thres))
    
    # plot(gene.expr.keep[,i], y)
    smoothScatter(x1.qc, y.qc, xlab = paste0("Gene expression for ", i), ylab = "chisq/df")
    # select the top5 influenctial points to avoid loss of power in the case of too many genes are dropped
    if(sum(cooksd > cooksd_thres) > 0){
      ig5 <- names(sort(cooksd, decreasing = TRUE)[1])
      # ig5 <- names(sort(cooksd, decreasing = TRUE)[1:pmin(5, sum(cooksd > cooksd_thres))])
    }else{
      ig5 <- NULL
    }
    remove.idx <- (names(cooksd) %in% ig5)
    keep.idx <- !(names(cooksd) %in% ig5)
    
    cat("\t Remove ", names(cooksd)[remove.idx], "\n")
    
    points(x1.qc[remove.idx], y.qc[remove.idx], col = "red")
    text(x = x1.qc, y = y.qc, labels = ifelse(remove.idx, names(cooksd), ""), col = "red")
    abline(lm0, col = "blue", lwd = 3, lty = 2)
    
    lm <- lm(y.qc[keep.idx] ~ x1.qc[keep.idx] + x2.qc[keep.idx], weights = df.qc[keep.idx]/2)
    abline(lm, col = "red", lwd = 3, lty = 2)
    chisq.pvals[i] <- exp(pt(summary(lm)$coefficients[2,3], df = summary(lm)$df[2], lower.tail = FALSE, log.p = TRUE))
    legend('top',legend = c(paste0('All genes: ', "p = ",  signif(chisq.pvals0[i], digits = 2)), 
                            paste0('Remove ', names(cooksd)[remove.idx], " with df = ", df.qc[remove.idx], ": p = ", signif(chisq.pvals[i], digits = 2))), 
           lwd = 3, lty = 2, col = c("blue", "red"), bty='n')
    
    j = j + 1
    
    if(is.null(ig5)){
      break  
    }
    # check again in a loop
    genes <- genes[keep.idx]
    y.qc <- y.qc[keep.idx]
    x1.qc <- x1.qc[keep.idx]
    x2.qc <- x2.qc[keep.idx]
    df.qc <- df.qc[keep.idx]
  }
}
dev.off()

chisq.pvals




q("no")






