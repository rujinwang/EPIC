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

trait = "blah_trait_blah"

rolypoly_out_dir <- file.path(curr_dir, "gtex_v8", "rolypoly_results", trait)
if(!dir.exists(rolypoly_out_dir)){
  dir.create(rolypoly_out_dir, recursive = TRUE)
}

load(file.path(input_dir, "gwasdata", "2_aligned", paste0(trait, "_aligned.rda")))
load(file.path(curr_dir, "snptogene", trait, paste0(trait, "_snp_to_gene.rda")))

annotation_dictionary <- read.table(file = file.path(curr_dir, "annotation_dictionary.txt"), stringsAsFactors = FALSE)
head(annotation_dictionary)
dim(annotation_dictionary)

loc <- read.table(file = file.path(curr_dir, "gtex_v8", paste0("gtex_v8.gene.noMHC.loc")), stringsAsFactors = FALSE)
head(loc)
dim(loc)

# sort the genes by genomic positions!!!
loc <- loc[order(loc$V2, loc$V3, loc$V4), ]
loc <- loc[which(!is.na(match(loc$V1, annotation_dictionary$V1))), ]
loc <- loc[which(!is.na(match(loc$V1, names(snp_to_gene)))), ]
dim(loc)

annotation_dictionary <- annotation_dictionary[which(!is.na(match(annotation_dictionary$V1, loc$V1))), ]
annotation_dictionary <- annotation_dictionary[which(!is.na(match(annotation_dictionary$V1, names(snp_to_gene)))), ]
dim(annotation_dictionary)

annotation_dictionary <- annotation_dictionary[match(loc$V1, annotation_dictionary$V1), ]
dim(annotation_dictionary)
head(annotation_dictionary)

length(snp_to_gene)
names(snp_to_gene)[1:5]

snp_to_gene <- snp_to_gene[match(annotation_dictionary$V1, names(snp_to_gene))]
length(snp_to_gene)
names(snp_to_gene)[1:5]

dim(gwas.align)
keep.snps <- unique(unlist(snp_to_gene, use.names = FALSE))
length(keep.snps)
head(keep.snps)
gwas.align <- gwas.align[match(keep.snps, gwas.align$rsid),]
dim(gwas.align)
head(gwas.align)

#################################################################################
# RolyPoly input: GWAS
gwas.align <- gwas.align[,c("chr", "pos", "rsid", "beta", "se", "MAF")]
dim(gwas.align)
head(gwas.align)

# Rename colnames of gwas summary stats
colnames(gwas.align) <- c("chrom", "pos", "rsid", "beta", "se", "maf")
dim(gwas.align)
head(gwas.align)

# RolyPoly input: Expression data
gene.expr = read.table(file = file.path(curr_dir, "gtex_v8", paste0("log2_tpm_median.txt")), header = TRUE, row.names = 1)
head(gene.expr); dim(gene.expr)
gene.expr.keep <- gene.expr[match(annotation_dictionary$V1, rownames(gene.expr)),]
head(gene.expr.keep); dim(gene.expr.keep)

# remove the average columns
gene.expr.keep <- gene.expr.keep[, !(colnames(gene.expr.keep) %in% "Average")]
head(gene.expr.keep); dim(gene.expr.keep)

# To compare across genes, we scaled and centered expression values across annotations.
gene.expr.keep.scaled = t(scale(t(gene.expr.keep)))
colnames(gene.expr.keep.scaled) = colnames(gene.expr.keep)
rownames(gene.expr.keep.scaled) = rownames(gene.expr.keep)
dim(gene.expr.keep.scaled); head(gene.expr.keep.scaled)

# According to the replies from Diego Calderon,
# we squared or took the absolute value of the gene expression values 
# (or normalized z-scores if I remember correctly) rather than 
# shifting them by the minimum
gene.expr.keep.scaled.abs = abs(gene.expr.keep.scaled)
dim(gene.expr.keep.scaled.abs); head(gene.expr.keep.scaled.abs)


# RolyPoly input: Gene annotation
# !!!!!!!!!!set stringAsFactors = FALSE, otherwise would get issues in spliting chromosomes
gene_annotation <- data.frame(chrom = annotation_dictionary$V2,
                              start = annotation_dictionary$V3,
                              end = annotation_dictionary$V4,
                              label = annotation_dictionary$V1,
                              stringsAsFactors = FALSE)
dim(gene_annotation)
head(gene_annotation)
tail(gene_annotation)

# copy from https://github.com/dcalderon/rolypoly/blob/master/R/data_io.R
# do the partition again!!!!!!!!!!!!!!!!!!!!!!
my_rolypoly_load_block_annotation <- function(rolypoly, block_annotation, genes = T) {
  message('adding block annotations')
  if (class(block_annotation) != 'data.frame') { stop('require data frame of annotations') }
  if (all(!(c('chrom', 'start', 'end', 'label') %in% colnames(block_annotation)))) {
    stop('require chrom, start, and end, cols')
  }
  if (any(duplicated(block_annotation$label))) { stop('remove duplicate block labels and rerun') }
  # sort and add a random partition label for bootstrap
  rolypoly$blocks <- block_annotation %>% dplyr::arrange(chrom, start) %>%
    dplyr::mutate(partition = cut(row_number(), breaks = 100, labels = F))
  return(rolypoly)
}

rp = list()
rp <- my_rolypoly_load_block_annotation(rp, block_annotation = gene_annotation)

ceiling(nrow(gene_annotation)/20)

if(trait == "LDL"){
  max.group = 418
  group.range = 1:max.group
} else if(trait == "TG"){
  max.group = 417
  group.range = 1:max.group
} else if(trait == "HDL"){
  max.group = 417
  group.range = 1:max.group
} else if(trait == "TC"){
  max.group = 418
  group.range = 1:max.group
} else if(trait == "SCZ"){
  max.group = 439 
  group.range = 1:max.group
} else if(trait == "BIP"){
  max.group = 439
  group.range = 1:max.group
} else if(trait == "SCZBIP"){
  # No group 184 for SCZBIP!!!
  max.group = 437
  group.range = c(1:183, 185:max.group)
} else if(trait == "T2Db"){
  max.group = 439 
  group.range = 1:max.group
}



rp.link <- NULL
for (group in group.range) {
  if(group %% 10 == 1){
    cat("Group = ", group, "\n") 
  }
  load(file = file.path(rolypoly_out_dir, "intermediate", paste0(trait, "_group", group, paste0("_rolypoly_link.rda"))))
  for (genei in names(rp.group$data)) {
    rp.group$data[[genei]]$partition <- rp$blocks$partition[rp$blocks$label == genei]
  }
  rp.link <- c(rp.link, rp.group$data)
}
length(rp.link)
# table(sapply(rp.link, function(z){z$partition}))

rp = list(data = rp.link)
nullgene = which(sapply(rp$data, function(z){length(z$partition)}) == 0)
rp$data = rp$data[-nullgene]


# Without any parallelization rolypoly could take a couple hours to run. 
# Much of the computation is spent linking SNPs to genes and reading in LD information. 
# If one provides the a rolypoly object to the rolypoly parameter of the rolypoly_roll 
# function call and a new object of block data, then only inference is performed. 
rp <- rolypoly_roll(
  rolypoly = rp, 
  gwas_data = gwas.align,
  block_annotation = gene_annotation,
  # block_data = gene.expr.keep.scaled.shift,
  block_data = gene.expr.keep.scaled.abs,
  bootstrap_iters = 100
)

# rp$bootstrap_results

save(rp, file = file.path(rolypoly_out_dir, paste0(trait, paste0("_rolypoly_inference.rda"))))


q("no")








