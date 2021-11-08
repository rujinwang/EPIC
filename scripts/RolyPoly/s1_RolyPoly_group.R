rm(list = ls())
curr_dir = "/pine/scr/r/u/rujin/EPIC"
input_dir = "/pine/scr/r/u/rujin/scRNAseqGWAS"

setwd(input_dir)
library(data.table)
library(rolypoly)
library(dplyr)
library(ggplot2)

trait = "blah_trait_blah"
group = blah_group_blah

# trait = "LDL"
# group = 1

rolypoly_out_dir <- file.path(curr_dir, "gtex_v8", "rolypoly_results", trait, "intermediate")
if(!dir.exists(rolypoly_out_dir)){
  dir.create(rolypoly_out_dir, recursive = TRUE)
}


load(file.path(input_dir, "gwasdata", "2_aligned", paste0(trait, "_aligned.rda")))
load(file.path(curr_dir, "snptogene", trait, paste0(trait, "_snp_to_gene.rda")))

annotation_dictionary <- read.table(file = file.path(curr_dir, "annotation_dictionary.txt"), stringsAsFactors = FALSE)
head(annotation_dictionary)
dim(annotation_dictionary)

loc <- read.table(file = file.path(curr_dir, "gtex_v8", paste0("gtex_v8.gene.loc")), stringsAsFactors = FALSE)
head(loc)
dim(loc)

# sort the genes by genomic positions!!!
loc <- loc[order(loc$V2, loc$V3, loc$V4), ]
loc <- loc[which(!is.na(match(loc$V1, annotation_dictionary$V1))), ]
loc <- loc[which(!is.na(match(loc$V1, names(snp_to_gene)))), ]
dim(loc)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!subset to speed up rolypoly linking procedure
start.idx <- (group - 1) * 20 + 1
end.idx <- pmin(group * 20, nrow(loc))

if(start.idx > end.idx){
  q("no")
}

loc <- loc[start.idx:end.idx, , drop = FALSE]


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

# RolyPoly input: Linkage disequilibrium (LD)
ld_path = file.path(input_dir, "rolypoly", "EUR_LD_FILTERED_NONAN_R")

rp.group <- rolypoly_roll(
  gwas_data = gwas.align,
  block_annotation = gene_annotation, 
  ld_folder = ld_path
)


save(rp.group, file = file.path(rolypoly_out_dir, paste0(trait, "_group", group, paste0("_rolypoly_link.rda"))))


q("no")








