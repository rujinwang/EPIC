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
library(nlshrink)

source(file = file.path(func_dir, "epic.R"))

type = "POET_QC"

trait = "blah_trait_blah"
chr = blah_chr_blah # 1:22
maf_thres = 0.01
sliding_step = 1
sliding_size = 10
Cmin = 1.5

load(file = file.path(curr_dir, "g1000_eur", paste0("genotype_chr", chr, ".rda")))
load(file = file.path(curr_dir, "g1000_eur", paste0("bim_chr", chr, ".rda")))
dim(X_ref)
dim(bim)
all(rownames(X_ref) == bim$id)

# Following Finucane et al., Nature Genetics, 2015, we removed SNPs within the major histocompatibility complex (MHC) region (Chr6: 25Mb- 34Mb)
loc <- read.table(file = file.path(curr_dir, "gtex_v8", paste0("gtex_v8.gene.noMHC.loc")), stringsAsFactors = FALSE)
head(loc)
dim(loc)
head(loc, 20)

# sort the genes by genomic positions!!!
loc <- loc[order(loc$V2, loc$V3, loc$V4), ]
loc <- loc[which(loc$V2==chr), ]
head(loc)
dim(loc)
tail(loc, 20)

out_dir <- file.path(curr_dir, "gtex_v8", "my_x_results", trait, paste0(type, "_", sliding_size, "_", Cmin, "_results"))
in_dir <- file.path(curr_dir, "gtex_v8", "my_x_results", trait, paste0(type, "_", sliding_size, "_", Cmin))
if(!dir.exists(out_dir)){
  dir.create(out_dir, recursive = TRUE)
}
if(!dir.exists(in_dir)){
  dir.create(in_dir, recursive = TRUE)
}

load(file = file.path(curr_dir, "trait_results", trait, paste0(trait, "_common_", maf_thres, "_chr", chr, "_chisq_spd_QC.rda")))
prunein_gene_pval <- prunein_gene_pval.chr[which(!is.na(match(names(prunein_gene_pval.chr), loc$V1)))]
length(prunein_gene_pval)

loc <- loc[which(!is.na(match(loc$V1, names(prunein_gene_pval)))), ]
dim(loc)

# reorder gene_pval list by position!!!!!!!!!!!!
prunein_gene_pval <- prunein_gene_pval[match(loc$V1, names(prunein_gene_pval))]
all(names(prunein_gene_pval)==loc$V1)

keep.snps <- unique(unlist(sapply(prunein_gene_pval, function(z){z$snps}), use.names = FALSE))
keep.ref.idx <- which(rownames(X_ref) %in% keep.snps)
X_ref <- X_ref[keep.ref.idx, ]
dim(X_ref)
bim <- bim[which(bim$id %in% keep.snps), ]
dim(bim)


C.list = unique(sort(c(Cmin, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)))
C.list = C.list[C.list >= Cmin]
C.list

gene.C <- rep(NA, length(prunein_gene_pval))
names(gene.C) <- names(prunein_gene_pval)

# calculate LD shrinkage estimators using sliding window!!!!!!!!!!

for (i in seq(1, (length(prunein_gene_pval) - sliding_size + 1), sliding_step)) {
  # snps in the sliding window
  gene.name <- names(prunein_gene_pval)[i]
  sw_snps <- unique(unlist(sapply(prunein_gene_pval[i:(i + sliding_size - 1)], function(z){z$snps}), use.names = FALSE))
  X_temp <- X_ref[which(rownames(X_ref) %in% sw_snps), ]
  
  cat("chr = ", chr, " gene ", i, ": #snps = ", nrow(X_temp), "\n")
  if(!file.exists(file.path(in_dir, paste0(trait, "_", maf_thres, "_LD_chr", chr, "_", gene.name, "_", type, ".txt"))) | 
     !file.exists(file.path(in_dir, paste0(trait, "_", maf_thres, "_LD_chr", chr, "_", gene.name, "_", type, ".rda")))){
    
    for (C in C.list) {
      cat("\t C = ", C, "\n")
      cov.POET <- POET(X_temp, K = 2, C = C, thres = 'soft', matrix = 'vad')
      corr.POET <- cov2cor(cov.POET$SigmaY)
      
      evs.POET <- eigen(corr.POET, only.values = TRUE)$values
      if(min(evs.POET) >= 0){
        break
      } else{
        if(C == 5){
          cat("ATTENTION!!!: CHECK THIS GENE")
        }
      }
    }
    gene.C[gene.name] = C
    
    write.table(names(prunein_gene_pval[i:(i + sliding_size - 1)]), file = file.path(in_dir, paste0(trait, "_", maf_thres, "_LD_chr", chr, "_", gene.name, "_", type, ".txt")), 
                quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
    
    save(cov.POET, corr.POET, gene.C, file = file.path(in_dir, paste0(trait, "_", maf_thres, "_LD_chr", chr, "_", gene.name, "_", type, ".rda")))
  }
}




q("no")





