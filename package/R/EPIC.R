#' @title Formatting GWAS summary statistics
#'
#' @description Large-scale GWAS summary statistics from different consortiums
#'  are in diverse formats, due to the nature of phenotypes being studied (case-control
#'  traits v.s. quantitative traits) and various analytical tools. We recommend formatting
#'  summary statistics in post-GWAS analysis. We perform quality control and allele
#'  alignment with respect to the external reference panel.
#'
#' @usage
#' format_gwas(gwas, bfile_path, plink_path)
#' @param gwas a data frame of GWAS summary statistics without formatting
#' @param bfile_path path to the \code{bim} file of reference panel
#' @param plink_path path to the PLINK toolkit
#'
#' @return A data frame of formatted GWAS summary statistics after quality control and allele alignment
#'
#' @examples
#' \dontrun{
#' bfile_path <- "/path/to/bfile"
#' plink_path <- "/path/to/PLINK"
#' gwas.formatted.Demo <- format_gwas(gwas.Demo, bfile_path, plink_path)
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import data.table glue
#' @export
format_gwas <- function(gwas, bfile_path, plink_path = NULL) {
  required.col <- c("chr", "pos", "rsid", "A1", "A2", "beta", "se", "P", "MAF", "N")
  if(!(all(required.col %in% colnames(gwas)))){
    col.index.lower <- which(tolower(colnames(gwas)) %in% c("chr", "pos", "rsid", "beta", "se"))
    col.index.upper <- which(toupper(colnames(gwas)) %in% c("A1", "A2", "P", "MAF", "N"))
    colnames(gwas)[col.index.lower] <- tolower(colnames(gwas))
    colnames(gwas)[col.index.upper] <- toupper(colnames(gwas))
  }
  if(!(all(required.col %in% colnames(gwas)))){
    stop("Required column of GWAS summary statistics is missing. ")
  }

  gwas$Zscore <- gwas$beta / gwas$se
  gwas$MAC <- round(2 * gwas$N * gwas$MAF)

  message("Chromosomes X, Y, and MT are filtered out...")
  keep.chr <- gwas$chr %in% seq_len(22)
  gwas <- gwas[keep.chr, ]

  message("Remove SNPs with MAF, MAC being NAs...")
  nonNA.snp <- which(!is.na(gwas$MAF) & !is.na(gwas$MAC))
  gwas <- gwas[nonNA.snp, ]
  dim(gwas)

  if(!file.exists(file.path(bfile_path, "plink.snplist"))){
    if(is.null(plink_path)){
      stop("Path to PLINK toolkit is missing. ")
    }
    setwd(bfile_path)
    cmd <- glue("{plink_path}/plink --bfile {bfile_path}/g1000_eur --maf 0.00000001 --geno 0 --write-snplist")
    message("Generate the list of SNP rs IDs in the external reference panel...")
    system(cmd)
  }
  message("Import the list of SNPs returned from PLINK...")
  snp.1000g <- fread(file.path(bfile_path, "plink.snplist"), header = FALSE)
  snp.1000g <- c(snp.1000g)
  snp.1000g <- snp.1000g[[1]]

  keep.snp <- gwas$rsid %in% snp.1000g
  gwas <- gwas[keep.snp, ]

  message("Remove mismatched SNPs with the reference panel...")
  ref <- fread(file = file.path(bfile_path, "g1000_eur.bim"), header = FALSE)
  ref.snps <- paste(ref$V1, ref$V4, ref$V5, ref$V6, sep = ":")
  snps <- paste(gwas$chr, gwas$pos, gwas$A1, gwas$A2, sep = ":")
  snps.flip <- paste(gwas$chr, gwas$pos, gwas$A2, gwas$A1, sep = ":")
  mismatch <- which(!(snps %in% ref.snps | snps.flip %in% ref.snps))
  if(length(mismatch)>0){
    gwas <- gwas[-mismatch,]
  }

  message("Align with the reference genome by flipping over the sign of beta and z-score...")
  snps.flip <- paste(gwas$chr, gwas$pos, gwas$A2, gwas$A1, sep = ":")

  flip.idx <- which(snps.flip %in% ref.snps)
  gwas.formatted <- gwas
  gwas.formatted$beta[flip.idx] <- -gwas$beta[flip.idx]
  gwas.formatted$Zscore[flip.idx] <- -gwas$Zscore[flip.idx]
  gwas.formatted$A1[flip.idx] <- gwas$A2[flip.idx]
  gwas.formatted$A2[flip.idx] <- gwas$A1[flip.idx]

  return(gwas.formatted)
}


#' @title Gene annotation
#'
#' @description We annotate each gene by defining a gene window with 10kb upstream of start site
#'  and 1.5kb downstream of end site, by default. One column of gene window length in kb is appended
#'  for future use in gene-wise LD pruning.
#'
#' @usage
#'  annotate_gene(upstream, downstream)
#' @param upstream numeric value of upsteam length from the start site of gene. Default is \code{10}.
#'  Unit is "kb".
#' @param downstream numeric value of downstream length from the end site of gene. Default is \code{1.5}.
#'  Unit is "kb".
#'
#' @return A data frame of gene annotations: gene name, chromosome, start site of annotation,
#'  end site of annotation, strand, and gene window size in kb.
#'
#' @examples
#' \dontrun{
#'  anno_dictionary <- annotate_gene(upstream = 10, downstream = 1.5)
#'  write.table(anno_dictionary, file = "annotation_dictionary.txt", quote = FALSE,
#'              sep = "\t", row.names = FALSE, col.names = FALSE)
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import utils
#' @export
annotate_gene <- function(upstream = 10, downstream = 1.5) {
  upstream = upstream * 1000
  downstream = downstream * 1000

  dictionary = read.table(file = system.file("extdata", "gene.noMHC.loc", package = "EPIC"), sep = "\t",
                          header = FALSE, stringsAsFactors = FALSE)

  anno_dictionary = dictionary
  for (i in seq_len(nrow(anno_dictionary))) {
    if(anno_dictionary$V5[i] == "+"){
      anno_dictionary$V3[i] = anno_dictionary$V3[i] - upstream
      anno_dictionary$V4[i] = anno_dictionary$V4[i] + downstream
    } else if(anno_dictionary$V5[i] == "-"){
      anno_dictionary$V3[i] = anno_dictionary$V3[i] - downstream
      anno_dictionary$V4[i] = anno_dictionary$V4[i] + upstream
    }
  }
  anno_dictionary = anno_dictionary[order(anno_dictionary$V2, anno_dictionary$V3),]

  window_kb = ceiling((anno_dictionary$V4 - anno_dictionary$V3)/1000)
  anno_dictionary = cbind(anno_dictionary, window_kb)

  return(anno_dictionary)
}


#' @title Mapping SNPs to genes
#'
#' @description By default, a gene window with 10kb upstream and 1.5kb downstream of each gene is defined
#'  and SNPs are assigned to genes based on genomic positional mapping.
#'
#' @usage
#'  map_snp_to_gene(gwas, anno_path = NULL)
#' @param gwas a data frame of formatted GWAS summary statistics
#' @param anno_path path to the user-specified annotation dictionary returned from \code{annotate_gene}.
#'  Default is 10kb upstream and 1.5kb downstream of each gene.
#'
#' @return A list of assigned SNPs for genes.
#'
#' @examples
#'  snp_to_gene.Demo <- map_snp_to_gene(gwas = gwas.formatted.Demo)
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import utils
#' @export
map_snp_to_gene <- function(gwas, anno_path = NULL) {
  if(is.null(anno_path)){
    anno_dictionary = read.table(file = system.file("extdata", "annotation_dictionary.txt", package = "EPIC"),
                                 sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  }else{
    anno_dictionary = read.table(file = file.path(anno_path),
                                 sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  }

  gwas.chrs <- unique(sort(gwas$chr))
  anno_dictionary <- anno_dictionary[anno_dictionary$V2 %in% gwas.chrs, ]
  snp_to_gene = vector('list', dim(anno_dictionary)[1])
  names(snp_to_gene) = anno_dictionary$V1
  for (i in seq_len(dim(anno_dictionary)[1])){
    if(i%%1000==1) {message("Assigning SNPs for gene ", i, "...")}
    gene.name = names(snp_to_gene)[i]
    chr = anno_dictionary$V2[i]
    start = anno_dictionary$V3[i]
    end = anno_dictionary$V4[i]
    idx = which(gwas$chr == chr & gwas$pos >= start & gwas$pos <= end)
    rsid = as.character(gwas[idx,"rsid"])
    snp_to_gene[[gene.name]] = rsid
  }

  zero.snps <- which(sapply(snp_to_gene, function(z){length(z)}) == 0)
  if(length(zero.snps)>0){
    snp_to_gene = snp_to_gene[-zero.snps]
  }
  return(snp_to_gene)
}


#' @title Divide SNPs into common and rare variants
#'
#' @description By default, EPIC use MAF=1% as the cutoff for common and rare variants.
#'  For rare variants, the upper bound of inclusion is controlled by MAF
#'  while the lower bound is determined by MAC (Minor Allele Count).
#'  Rare variants with MAC less than 20 are removed from analysis.
#'
#' @usage
#'  divide_common_rare(gwas, snp_to_gene, maf_thres = 0.01, mac_lb = 20)
#' @param gwas a data frame of formatted GWAS summary statistics
#' @param snp_to_gene A list of assigned SNPs for genes, returned from \code{map_snp_to_gene}.
#' @param maf_thres the MAF cutoff for common and rare variants. Default is \code{0.01}.
#' @param mac_lb the minor allele count (MAC) lower bound. A gene is excluded from the analysis if its MAC is smaller
#'  than the MAC lower bound. Default is \code{20}.
#'
#' @return A list with components
#'  \item{common_snp_to_gene}{A list of common variants for genes}
#'  \item{rare_snp_to_gene}{A list of rare variants for genes}
#'
#' @examples
#' \dontrun{
#' snp_division_obj.Demo <- divide_common_rare(gwas = gwas.formatted.Demo,
#'                                             snp_to_gene = snp_to_gene.Demo)
#' common_snp_to_gene.Demo <- snp_division_obj.Demo$common_snp_to_gene
#' rare_snp_to_gene.Demo <- snp_division_obj.Demo$rare_snp_to_gene
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @export
divide_common_rare <- function(gwas, snp_to_gene, maf_thres = 0.01, mac_lb = 20) {
  if(maf_thres > 1 | maf_thres < 0 | mac_lb <= 0){
    stop("Invalid MAF threshold or MAC lower bound! ")
  }

  rare.snps <- as.character(gwas$rsid[gwas$MAF <= maf_thres])
  gwas.rare <- gwas[which(gwas$MAF <= maf_thres),]

  # Map rare and common SNPs to genes
  rare_snp_to_gene = vector('list', length(snp_to_gene))
  names(rare_snp_to_gene) = names(snp_to_gene)
  for (i in seq_len(length(names(rare_snp_to_gene)))){
    if(i%%200==1) {message("Assigning rare variants for gene ", i)}
    rare_snp_to_gene[[i]] = snp_to_gene[[i]][which(snp_to_gene[[i]] %in% rare.snps)]
  }
  message("\n", "Assigning common variants for genes...")
  common_snp_to_gene <- mapply(function(x,y) {setdiff(x, y)}, x = snp_to_gene, y = rare_snp_to_gene)

  zero.rare.genes <- which(sapply(rare_snp_to_gene, function(z){length(z)}) == 0)
  if(length(zero.rare.genes)>0){
    rare_snp_to_gene <- rare_snp_to_gene[-zero.rare.genes]
  }

  zero.common.genes <- which(sapply(common_snp_to_gene, function(z){length(z)}) == 0)
  if(length(zero.common.genes)>0){
    common_snp_to_gene <- common_snp_to_gene[-zero.common.genes]
  }

  message("\n", "Remove genes with low MAC for rare variants... ")
  # gene.MAC <- sapply(rare_snp_to_gene, function(z){sum(gwas.rare$MAC[match(z, gwas.rare$rsid)])})
  # gene.low.MAC <- which(gene.MAC < mac_lb)
  # if(length(gene.low.MAC)>0){
  #   rare_snp_to_gene.test <- rare_snp_to_gene[-gene.low.MAC]
  # }
  snps.low.MAC <- as.character(gwas.rare$rsid[gwas.rare$MAC < mac_lb])
  gene.with.low.MAC <- sapply(rare_snp_to_gene, function(z){any(z %in% snps.low.MAC)})
  gene.low.MAC <- rep(FALSE, length(rare_snp_to_gene))
  gene.MAC <- sapply(rare_snp_to_gene[gene.with.low.MAC], function(z){sum(gwas.rare$MAC[match(z, gwas.rare$rsid)])})
  index.gene.MAC <- which(gene.MAC < mac_lb)
  gene.low.MAC[gene.with.low.MAC][index.gene.MAC] = TRUE
  if(sum(gene.low.MAC)>0){
    rare_snp_to_gene <- rare_snp_to_gene[!gene.low.MAC]
  }

  return(list(common_snp_to_gene = common_snp_to_gene, rare_snp_to_gene = rare_snp_to_gene))
}


#' @title Perform LD pruning on common variants and import prune-in variants
#'
#' @description EPIC imports prune-in variants following LD pruning pipeline for
#' common variants using PLINK.
#'
#' @usage
#'  LD_pruning(common_snp_to_gene, prunein_dir)
#' @param common_snp_to_gene A list of assigned common SNPs for genes, returned from \code{divide_common_rare}.
#' @param prunein_dir the directory of LD pruning intermediate files
#'
#' @return A list of assigned prune-in common SNPs for genes.
#'
#' @examples
#' \dontrun{
#' out_dir = "/path/to/plink_output"
#' MAF <- 0.01
#' prunein_dir <- file.path(out_dir, "prunein", MAF)
#' prunein_snp_to_gene.Demo <- LD_pruning(common_snp_to_gene.Demo, prunein_dir)
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import data.table
#' @export
LD_pruning <- function(common_snp_to_gene, prunein_dir) {
  prunein_snp_to_gene = vector('list', length(common_snp_to_gene))
  names(prunein_snp_to_gene) = names(common_snp_to_gene)

  for (i in seq_len(length(names(common_snp_to_gene)))){
    if(i%%100==1) {message("\n", "Assigning prune-in common variants for gene", i)}
    gene = names(common_snp_to_gene)[i]
    if(file.exists(file.path(prunein_dir, paste0(gene, ".prune.in")))){
      prunein.file = fread(file.path(prunein_dir, paste0(gene, ".prune.in")), header = FALSE)
      prunein_snp_to_gene[[i]] = common_snp_to_gene[[i]][(common_snp_to_gene[[i]] %in% prunein.file$V1)]
    }
  }

  zero.prunein.genes <- which(sapply(prunein_snp_to_gene, function(z){length(z)}) == 0)
  if(length(zero.prunein.genes)>0){
    prunein_snp_to_gene <- prunein_snp_to_gene[-zero.prunein.genes]
  }
  return(prunein_snp_to_gene)
}





#' @title Perform LD pruning on common variants and import prune-in variants
#'
#' @description EPIC imports prune-in variants following LD pruning pipeline for
#' common variants using PLINK.
#'
#' @usage
#'  second_pruning(prunein_snp_to_gene, plink_path, bfile_path, chr, prunein_dir)
#' @param prunein_snp_to_gene A list of assigned common SNPs for genes, returned from \code{LD_pruning}.
#' @param bfile_path path to the \code{bim} file of reference panel
#' @param plink_path path to the PLINK toolkit
#' @param chr numeric value of chromosome number
#' @param prunein_dir the directory of LD pruning intermediate files
#'
#' @return A list of assigned prune-in common SNPs for genes.
#'
#' @examples
#' \dontrun{
#' out_dir = "/path/to/plink_output"
#' trait = "LDL"
#' MAF = 0.01
#' prunein_dir <- file.path(out_dir, paste0("prunein2", "_", trait), MAF)
#' bfile_path <- "/path/to/bfile"
#' plink_path <- "/path/to/PLINK"
#' chr = 2
#' prunein_qc_snp_to_gene.chr <- second_pruning(prunein_snp_to_gene.Demo,
#'                                              plink_path, bfile_path, chr, prunein_dir)
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import data.table glue utils
#' @export
second_pruning <- function(prunein_snp_to_gene, plink_path, bfile_path, chr, prunein_dir) {
  anno_dictionary <- read.table(file = system.file("extdata", "annotation_dictionary.txt", package = "EPIC"),
                                sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gene.chr <- as.character(anno_dictionary$V1)[anno_dictionary$V2==chr]
  prunein_qc_snp_to_gene <- prunein_qc_snp_to_gene[which(names(prunein_qc_snp_to_gene) %in% gene.chr)]
  prunein.size <- sapply(prunein_qc_snp_to_gene, function(z){length(z)}, USE.NAMES = FALSE)

  for (i in seq_len(length(prunein_qc_snp_to_gene))) {
    gene.name = names(prunein_qc_snp_to_gene)[i]

    # gene-window size
    gwindow <- anno_dictionary$V6[which(as.character(anno_dictionary$V1)==gene.name)]
    gfrom <- anno_dictionary$V3[which(as.character(anno_dictionary$V1)==gene.name)]
    gto <- anno_dictionary$V4[which(as.character(anno_dictionary$V1)==gene.name)]
    gchr <- anno_dictionary$V2[which(as.character(anno_dictionary$V1)==gene.name)]

    ld_thres <- 0.8
    while(prunein.size[i] > 200){
      message("Second-round LD pruning on ", i, ": ", gene.name)

      if(prunein.size[i] > 2000){
        ld_thres <- 0.05
      } else{
        ld_thres <- ld_thres - 0.1
      }

      if(ld_thres < 0.1){
        ld_thres <- 0.05
      }

      cat("LD pruning at ", ld_thres, "\n")

      if(!dir.exists(file.path(prunein_dir, ld_thres))){
        dir.create(file.path(prunein_dir, ld_thres), recursive = TRUE)
      }
      out <- file.path(prunein_dir, ld_thres)


      if(!file.exists(file.path(out, paste0(gene.name, ".prune.in")))){
        temp.file <- tempfile(gene.name, fileext = ".txt")
        write.table(prunein_qc_snp_to_gene[[gene.name]], file = temp.file,
                    quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)

        cmd <- glue("{plink_path} --bfile {bfile_path} --maf {maf_thres} --from-bp {gfrom} --to-bp {gto} --chr {gchr} --indep-pairwise {gwindow}'kb' 1 {ld_thres} --extract {temp.file} --out {out}/{gene.name}")
        system(cmd)
      }

      prunein.file = fread(file.path(out, paste0(gene.name, ".prune.in")), header = FALSE)
      snps <- prunein_qc_snp_to_gene[[gene.name]]
      prunein_snps <- prunein.file$V1

      match.idx <- which(!is.na(match(snps, prunein_snps)))
      prunein_qc_snp_to_gene[[gene.name]] <- prunein_qc_snp_to_gene[[gene.name]][match.idx]
      prunein.size[i] <- length(prunein_qc_snp_to_gene[[gene.name]])

      if(ld_thres < 0.1){
        break
      }
    }
  }
  return(prunein_qc_snp_to_gene)
}





#' @title GWAS summary statistics importing
#'
#' @description Import GWAS summary statisics for gene
#'
#' @usage
#'  read_in(gwas, snp_to_gene)
#' @param gwas a data frame of formatted GWAS summary statistics
#' @param snp_to_gene A list of assigned SNPs for genes
#'
#' @return A list of GWAS summary statistics for genes.
#'
#' @examples
#' prunein_gene_pval.Demo <- read_in(gwas = gwas.formatted.Demo,
#'                                   snp_to_gene = prunein_qc_snp_to_gene.Demo)
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import data.table
#' @export
read_in <- function(gwas, snp_to_gene) {
  mapped.snps <- unique(unlist(snp_to_gene, use.names = FALSE))
  gwas <- gwas[match(mapped.snps, gwas$rsid),]
  gwas <- as.data.table(gwas)

  gene_pval = vector("list", length = length(snp_to_gene))
  names(gene_pval) = names(snp_to_gene)

  for (i in seq_len(length(gene_pval))) {
    gene.name = names(gene_pval)[i]
    if(i%%200==1) {message("Read in GWAS summary statistics for gene ", i, ": ", gene.name)}
    gene_pval[[gene.name]] = vector("list", 5)
    names(gene_pval[[gene.name]]) = c("snps", "betahat", "sehat", "zhat", "gwasP")
    idx = match(snp_to_gene[[i]], gwas$rsid)
    gene_pval[[gene.name]][["snps"]] = as.character(gwas$rsid[idx])
    gene_pval[[gene.name]][["betahat"]] = gwas$beta[idx]
    gene_pval[[gene.name]][["sehat"]] = gwas$se[idx]
    gene_pval[[gene.name]][["zhat"]] = gwas$Zscore[idx]
    gene_pval[[gene.name]][["gwasP"]] = gwas$P[idx]
  }
  return(gene_pval)
}


#' @title Transcriptomic data processing for scRNA-seq
#'
#' @description For scRNA-seq data, we follow the Seurat pipeline to perform gene- and cell-wise
#' quality controls and focus on the top 8000 highly variable genes.
#' Cell-type-specific RPKMs are calculated by combining reads or UMI counts from
#' all cells pertaining to a specific cell type, followed by log2 transformation
#' with an added pseudo-count.
#'
#' @usage
#'  process_scRNA(SeuratObject, meta_ct_col)
#' @param SeuratObject a Seurat object from a gene expression matrix, returned from \code{CreateSeuratObject}.
#' @param meta_ct_col column name of cell types in the meta data matrix.
#'
#' @return A list with components
#'  \item{scRNA.rpkm}{A data frame of cell-type-specific RPKMs with log2 transformation}
#'  \item{scRNA.loc}{A data frame of feature dictionary for top 8000 highly variable genes}
#'
#' @examples
#' \dontrun{
#'  scRNA.object <- process_scRNA(SeuratObject, meta_ct_col = "cell.type")
#'  scRNA.rpkm <- scRNA.object$scRNA.rpkm
#'  scRNA.loc <- scRNA.object$scRNA.loc
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import Seurat biomaRt utils
#' @importFrom  BiocGenerics t
#' @export
process_scRNA <- function(SeuratObject, meta_ct_col) {
  gene.loc = read.table(file = system.file("extdata", "gene.noMHC.loc", package = "EPIC"), sep = "\t",
                        header = FALSE, stringsAsFactors = FALSE)

  grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice")
  grch37 = useDataset("hsapiens_gene_ensembl", mart=grch37)

  HGNC_to_ENSE=getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', "strand"),
                     filters = 'hgnc_symbol',
                     values = rownames(SeuratObject),
                     mart = grch37)

  dictionary = as.data.frame(HGNC_to_ENSE)
  dictionary <- dictionary[dictionary$ensembl_gene_id %in% gene.loc$V1, ]

  dictionary$chromosome_name <- as.numeric(dictionary$chromosome_name)
  o <- order(dictionary$chromosome_name, dictionary$start_position, dictionary$end_position)
  dictionary <- dictionary[o,]

  rawCount <- GetAssayData(object = SeuratObject, assay = "RNA", slot = "counts")
  keep.idx <- which(!is.na(match(rownames(rawCount), dictionary$hgnc_symbol)))
  rawCount <- rawCount[keep.idx, ]
  rawCount <- rawCount[match(dictionary$hgnc_symbol, rownames(rawCount)), ]

  # RPKM
  gene.length <- (dictionary$end_position - dictionary$start_position)/1000
  library.size <- apply(rawCount, 2, sum)
  rpkm <- (10^6)*t(t(rawCount)/library.size)/gene.length
  cts <- sort(unique(SeuratObject@meta.data[[meta_ct_col]]))

  log2_rpkm_ct <- matrix(NA, nrow = nrow(rpkm), ncol = length(cts))
  rownames(log2_rpkm_ct) <- rownames(rpkm)
  colnames(log2_rpkm_ct) <- cts
  for (ct in cts) {
    cell.idx = rownames(SeuratObject@meta.data)[which(SeuratObject@meta.data[[meta_ct_col]] == ct)]
    log2_rpkm_ct[,ct] <- log2(1 + apply(rpkm[, cell.idx, drop = FALSE], 1, mean))
  }

  SeuratObject <- subset(SeuratObject, features = rownames(rpkm))

  ngene = 8000
  var.features <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = ngene)

  common.genes.ngene <- var.features@assays$RNA@var.features
  dictionary.ngene <- dictionary[dictionary$hgnc_symbol %in% common.genes.ngene, ]
  log2_rpkm_ct.ngene <- log2_rpkm_ct[rownames(log2_rpkm_ct) %in% common.genes.ngene, ]
  colnames(log2_rpkm_ct.ngene) <- gsub(" ", "_", colnames(log2_rpkm_ct.ngene))

  log2_rpkm_ct.ngene = cbind(log2_rpkm_ct.ngene, apply(log2_rpkm_ct.ngene, 1, mean))
  colnames(log2_rpkm_ct.ngene)[ncol(log2_rpkm_ct.ngene)] = "Average"
  rownames(log2_rpkm_ct.ngene) <- dictionary$ensembl_gene_id[match(rownames(log2_rpkm_ct.ngene), dictionary$hgnc_symbol)]

  log2_rpkm_ct.ngene = as.data.frame(log2_rpkm_ct.ngene)

  dictionary.ngene <- dictionary.ngene[,-1]
  colnames(dictionary.ngene) <- c("V1","V2","V3","V4","V5")

  return(list(scRNA.rpkm = log2_rpkm_ct.ngene, scRNA.loc = dictionary.ngene))
}




#' @title POET estimator and sliding window approach
#'
#' @description We adopt the POET estimator to obtain a well-conditioned covariance matrix
#'  via sparse shrinkage. We implement a sliding-window approach to only compute
#'  correlations for pairs of genes within a certain distance from each.
#'
#' @usage
#'  calculate_POET_sw(genotype, gene_pval, gene.loc, chr, type, inter_dir,
#'                    X_super = NULL, sliding_size = 10, Cmin = 1.5, sliding_step = 1)
#'
#' @param genotype a genotype matrix
#' @param gene_pval a list of GWAS summary statistics for genes
#' @param gene.loc a data frame of gene dictionary, only including genes of interest in each transcriptomic profile
#'  (e.g., top 8000 highly variable genes from scRNA-seq or genes with great specificity score from GTEx v8).
#'  This can be returned from \code{process_scRNA}.
#' @param chr numeric value of chromosome number
#' @param type type of variants included in POET estimation. This should be either 'prunein' or 'joint'.
#' @param inter_dir directory of intermediate output files
#' @param X_super a matrix of aggregated pseudo-SNP. Default is \code{NULL}.
#' @param sliding_size size of sliding window. Default is \code{10}.
#' @param Cmin a user-specified minimum constant of thresholding funciton, determining the POET shrinkage.
#'  Default is \code{1.5}.
#' @param sliding_step the moving step of sliding-window. Default is \code{1}.
#'
#' @return Intermediate stored files including
#'  \item{.rda}{POET shrinkage correlation estimator}
#'  \item{.txt}{A list of genes residing in the same sliding window}
#'
#' @examples
#' \dontrun{
#' inter_dir <- file.path("/path/to/intermediates", trait)
#' calculate_POET_sw(genotype = X_ref.Demo, gene_pval = prunein_gene_pval.Demo,
#'                   gene.loc = gtex.loc,
#'                   chr = chr, type = "POET", inter_dir = inter_dir)
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats utils
#' @export
calculate_POET_sw <- function(genotype, gene_pval, gene.loc, chr, type, inter_dir, X_super = NULL, sliding_size = 10, Cmin = 1.5, sliding_step = 1) {
  if(!type %in% c("POET", "iPOET")){
    stop("Type should be either POET or iPOET! ")
  }
  if(Cmin <=0){
    stop("Cmin should be a positive numeric value! ")
  }
  gene.loc <- gene.loc[order(gene.loc$V2, gene.loc$V3, gene.loc$V4), ]
  gene.loc <- gene.loc[which(gene.loc$V2==chr), ]

  inter_dir <- file.path(inter_dir, type)
  if(!dir.exists(inter_dir)){
    dir.create(inter_dir, recursive = TRUE)
  }
  gene.loc <- gene.loc[which(!is.na(match(gene.loc$V1, names(gene_pval)))), ]
  gene_pval <- gene_pval[match(gene.loc$V1, names(gene_pval))]
  keep.snps <- unique(unlist(sapply(gene_pval, function(z){z$snps}), use.names = FALSE))
  keep.ref.idx <- which(rownames(genotype) %in% keep.snps)
  genotype <- genotype[keep.ref.idx, ]

  C.list = unique(sort(c(Cmin, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)))
  C.list = C.list[C.list >= Cmin]

  for (i in seq(1, (length(gene_pval) - sliding_size + 1), sliding_step)) {
    gene.name <- names(gene_pval)[i]
    sw_snps <- unique(unlist(sapply(gene_pval[i:(i + sliding_size - 1)], function(z){z$snps}), use.names = FALSE))
    if(type == "iPOET" & !is.null(X_super)){
      X_temp <- genotype[which(rownames(genotype) %in% sw_snps), ]
      super_genes <- names(gene_pval)[i:(i + sliding_size - 1)]
      X_rare <- X_super[which(rownames(X_super) %in% super_genes), , drop = FALSE]
      if(nrow(X_rare)==0){
        X_temp <- X_temp
      } else{
        X_temp <- rbind(X_temp, X_rare)
      }
    }else{
      X_temp <- genotype[which(rownames(genotype) %in% sw_snps), ]
    }
    message("chr = ", chr, " gene ", i, ": #snps = ", nrow(X_temp))
    if(!file.exists(file.path(inter_dir, paste0("LD_chr", chr, "_", gene.name, "_", type, ".txt"))) |
       !file.exists(file.path(inter_dir, paste0("LD_chr", chr, "_", gene.name, "_", type, ".rda")))){
      for (C in C.list) {
        message("\t C = ", C)
        cov.POET <- POET(X_temp, K = 2, C = C, thres = 'soft', matrix = 'vad')
        corr.POET <- cov2cor(cov.POET$SigmaY)
        evs.POET <- eigen(corr.POET, only.values = TRUE)$values
        if(min(evs.POET) >= 0){
          break
        } else{
          if(C == 5){
            message("ATTENTION!!!: CHECK THIS GENE")
          }
        }
      }
      write.table(names(gene_pval[i:(i + sliding_size - 1)]), file = file.path(inter_dir, paste0("LD_chr", chr, "_", gene.name, "_", type, ".txt")),
                  quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
      save(corr.POET, file = file.path(inter_dir, paste0("LD_chr", chr, "_", gene.name, "_", type, ".rda")))
    }
  }
}



POET <- function(Y, K, C, thres = 'soft', matrix = 'vad') {
  p <- nrow(Y)
  n <- ncol(Y)
  Y <- Y - t(t(apply(t(Y), 2, mean)))%*%matrix(1, 1, n)

  if(K > 0){
    V <- eigen(t(Y)%*%Y)$vectors
    temp = try(V > 0, silent=TRUE)
    if(is.character(temp)){
      V <- Re(V)
    }

    V <- as.matrix(V)
    Dd <- eigen(t(Y)%*%Y)$values
    Dd <- as.vector(Dd)
    W <- sort(diag(Dd),index.return=TRUE)$x
    W <- as.matrix(W)
    Id <- sort(diag(Dd),index.return=TRUE)$ix
    Id <- as.matrix(Id)
    F <- sqrt(n)*V[,seq_len(K)]
    LamPCA <- Y%*%F/n
    uhat <- Y-LamPCA%*%t(F)
    Lowrank <- LamPCA%*%t(LamPCA)
    rate <- 1/sqrt(p)+sqrt((log(p))/n)
  } else {
    uhat <- Y
    rate <- sqrt((log(p))/n)
    Lowrank <- matrix(0,p,p)
  }

  SuPCA <- uhat%*%t(uhat)/n
  SuDiag <- diag(diag(SuPCA))
  if(matrix == 'cor'){
    R <- solve(SuDiag^(1/2))%*%SuPCA%*%solve(SuDiag^(1/2))
  }
  if(matrix == 'vad') {
    R <- SuPCA
  }

  uu <- array(0,dim=c(p,p,n))
  roottheta <- array(0,dim=c(p,p))
  lambda <- array(0,dim=c(p,p))
  for (i in seq_len(p)){
    for (j in seq_len(i)){
      uu[i,j,] <- uhat[i,]*uhat[j,]
      roottheta[i,j] <- sd(uu[i,j,])
      lambda[i,j] <- roottheta[i,j]*rate*C
      lambda[j,i] <- lambda[i,j]
    }
  }

  Rthresh=matrix(0,p,p)

  if(thres == 'soft'){
    for (i in seq_len(p)){
      for (j in seq_len(i)){
        if (abs(R[i,j])<lambda[i,j] && j<i){
          Rthresh[i,j] <- 0
        }else{
          if (j==i){
            Rthresh[i,j] <- R[i,j]
          }else{
            Rthresh[i,j] <- sign(R[i,j])*(abs(R[i,j])-lambda[i,j])
          }
        }
        Rthresh[j,i] <- Rthresh[i,j];
      }
    }
  }

  SigmaU <- matrix(0,p,p)
  if (matrix == 'cor'){
    SigmaU <- SuDiag^(1/2)%*%Rthresh*SuDiag^(1/2)
  }
  if (matrix == 'vad')     {
    SigmaU <- Rthresh
  }
  SigmaY <- SigmaU+Lowrank

  return(list(SigmaU=SigmaU,SigmaY=SigmaY,factors=t(F),loadings=LamPCA))
}





if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("corr.POET"))
}

#' @title Chi-square gene-level association testing and gene-gene correlation
#'
#' @description We perform gene-level chi-square association testing and calculate
#'  gene-gene correlations.
#'
#' @usage
#'  get_gene_chisq(gene_pval, gene.loc, type = "POET", inter_dir,
#'                 max.dist = 5, sliding_size = 10, sliding_step = 1)
#'
#' @param gene_pval a list of GWAS summary statistics for genes
#' @param gene.loc a data frame of gene dictionary, only including genes of interest in each transcriptomic profile
#'  (e.g., top 8000 highly variable genes from scRNA-seq or genes with great specificity score from GTEx v8).
#'  This can be returned from \code{process_scRNA}.
#' @param type type of variants included in POET estimation. This should be either 'prunein' or 'joint'.
#' @param inter_dir directory of intermediate output files
#' @param max.dist the maximum distance within which pair-wise gene-gene correlation is computed. Default is \code{5}. Unit is "Mb".
#' @param sliding_size size of sliding window. Default is \code{10}.
#' @param sliding_step the moving step of sliding-window. Default is \code{1}.
#'
#' @return A list with components
#'  \item{gene_POET_pval}{A list of gene-level testing}
#'  \item{gene_corr.mat}{A matrix of gene-gene correlations}
#'
#' @examples
#' \dontrun{
#' pruneinObject.Demo <- get_gene_chisq(gene_pval = prunein_gene_pval.Demo,
#'                                      gene.loc = gtex.loc,
#'                                      type = "POET", inter_dir = inter_dir)
#' prunein_POET_gene_pval.Demo <- pruneinObject.Demo$gene_POET_pval
#' prunein_gene_corr.mat.Demo <- pruneinObject.Demo$gene_corr.mat
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats utils
#' @export
get_gene_chisq <- function(gene_pval, gene.loc, type = "POET", inter_dir, max.dist = 5, sliding_size = 10, sliding_step = 1) {
  if(!type %in% c("POET", "iPOET")){
    stop("Type should be either POET or iPOET! ")
  }
  max.dist.bp = max.dist * 1000000
  inter_dir <- file.path(inter_dir, type)

  gene.loc <- gene.loc[order(gene.loc$V2, gene.loc$V3, gene.loc$V4), ]
  gene_pval <- gene_pval[which(!is.na(match(names(gene_pval), gene.loc$V1)))]
  gene.loc <- gene.loc[which(!is.na(match(gene.loc$V1, names(gene_pval)))), ]
  gene_pval <- gene_pval[match(gene.loc$V1, names(gene_pval))]

  gene_corr = vector("list", length = length(gene_pval))
  names(gene_corr) = names(gene_pval)

  gene_POET_pval <- vector('list', length = length(gene_pval))
  names(gene_POET_pval) <- names(gene_pval)

  for(chr in intersect(seq_len(22), unique(gene.loc$V2))){
    nsnps <- sum(gene.loc$V2==chr)
    genes.list <- gene.loc$V1[which(gene.loc$V2==chr)]
    gene.loc.temp <- gene.loc[which(gene.loc$V2==chr), ]
    for (i in seq(1, (nsnps - sliding_size + 1), sliding_step)) {
      load(file = file.path(inter_dir, paste0("LD_chr", chr, "_", genes.list[i], "_", type, ".rda")))
      sw.gene <- c(read.table(file.path(inter_dir, paste0("LD_chr", chr, "_", genes.list[i], "_", type, ".txt")),
                              stringsAsFactors = FALSE)$V1)
      gene.name = genes.list[i]

      message(i, ": ", gene.name, " chr = ", chr)

      if(is.null(gene_pval[[gene.name]]$rare_snps)){
        gene_POET_pval[[gene.name]]$snps = gene_pval[[gene.name]]$snps
        snps = gene_POET_pval[[gene.name]]$snps
        gene_POET_pval[[gene.name]]$zhat = gene_pval[[gene.name]]$zhat
        gene_POET_pval[[gene.name]]$gwasP = gene_pval[[gene.name]]$gwasP
      }else{
        gene_POET_pval[[gene.name]]$snps = c(gene_pval[[gene.name]]$snps, gene.name)
        snps = gene_POET_pval[[gene.name]]$snps
        gene_POET_pval[[gene.name]]$zhat = c(gene_pval[[gene.name]]$zhat, gene_pval[[gene.name]]$rare_T)
        gene_POET_pval[[gene.name]]$gwasP = c(gene_pval[[gene.name]]$gwasP, gene_pval[[gene.name]]$rare_pval)
      }

      snps = gene_POET_pval[[gene.name]]$snps
      gene_POET_pval[[gene.name]]$LD <- corr.POET[snps, snps, drop = FALSE]
      LD_inverse = compute_ld_inverse(gene_POET_pval[[gene.name]][["LD"]], tol = 1e-5)
      gene_POET_pval[[gene.name]]$inv_LD = LD_inverse$inv_mat
      gene_POET_pval[[gene.name]]$eigenvalue = LD_inverse$ev
      gene_POET_pval[[gene.name]]$df = LD_inverse$df
      z_mat = as.matrix(gene_POET_pval[[gene.name]][["zhat"]])
      gene_POET_pval[[gene.name]]$chisq = as.numeric(t(z_mat) %*% gene_POET_pval[[gene.name]][["inv_LD"]] %*% z_mat)
      gene_POET_pval[[gene.name]]$pval = exp(pchisq(gene_POET_pval[[gene.name]][["chisq"]], df = gene_POET_pval[[gene.name]][["df"]], log.p = TRUE, lower.tail = FALSE))

      # gene-gene correlation
      gene1 = gene.name
      gene1.snps <- gene_POET_pval[[gene1]][["snps"]]

      st <- gene.loc$V3[as.character(gene.loc$V1)==gene1]
      ed <- gene.loc$V4[as.character(gene.loc$V1)==gene1]

      temp <- gene.loc.temp[-seq_len(i), , drop = FALSE]
      close.gene <- as.character(temp$V1[temp$V2 == chr & temp$V3 < ed + max.dist.bp & temp$V4 > st - max.dist.bp]) # 5MB
      close.gene <- close.gene[close.gene %in% sw.gene]
      close.gene <- close.gene[close.gene %in% names(gene_POET_pval)]

      if(length(close.gene)>0){
        gene_corr[[gene1]] = vector("list", length = length(close.gene))
        names(gene_corr[[gene1]]) = close.gene

        if(length(close.gene)>0){
          for(gene2 in close.gene){
            if(is.null(gene_pval[[gene2]]$rare_snps)){
              gene2.snps <- gene_pval[[gene2]][["snps"]]
            } else{
              gene2.snps <- c(gene_pval[[gene2]][["snps"]], gene2)
            }

            R <- corr.POET[c(gene1.snps, gene2.snps), c(gene1.snps, gene2.snps), drop = FALSE]

            p <- length(gene1.snps)
            q <- length(gene2.snps)
            gene1.invR.sqrt = compute_ld_inverse_sqrt(R[seq_len(p), seq_len(p), drop = FALSE], tol = 1e-5)
            gene2.invR.sqrt = compute_ld_inverse_sqrt(R[(p+1):(p+q), (p+1):(p+q), drop = FALSE], tol = 1e-5)

            three.parts = gene1.invR.sqrt %*% R[seq_len(p), (p+1):(p+q)] %*% gene2.invR.sqrt
            R.star <- diag(p+q)
            R.star[seq_len(p), (p+1):ncol(R.star)] <- three.parts
            R.star[(p+1):ncol(R.star), seq_len(p)] <- t(three.parts)

            R.star.chol <- cholesky_decomposition(R.star)
            cov.chisq <- compute_covariance_cholesky(Lmat = R.star.chol, p = p, q = q)
            corr.chisq <- cov.chisq / (sqrt(2*p) * sqrt(2*q))

            gene_corr[[gene1]][[gene2]]$R <- R
            gene_corr[[gene1]][[gene2]]$cov.chisq <- cov.chisq
            gene_corr[[gene1]][[gene2]]$corr.chisq <- corr.chisq

            # message("\t corr with ", gene2, "=", corr.chisq)
          }
        }
      }
    }

    lastn <- length(genes.list) - sliding_size + 2
    load(file = file.path(inter_dir, paste0("LD_chr", chr, "_", genes.list[(lastn - 1)], "_", type, ".rda")))
    for (i in seq(lastn, length(genes.list), sliding_step)) {
      gene.name = genes.list[i]
      message(i, ": ", gene.name, " chr = ", chr)

      if(is.null(gene_pval[[gene.name]]$rare_snps)){
        gene_POET_pval[[gene.name]]$snps = gene_pval[[gene.name]]$snps
        snps = gene_POET_pval[[gene.name]]$snps
        gene_POET_pval[[gene.name]]$zhat = gene_pval[[gene.name]]$zhat
        gene_POET_pval[[gene.name]]$gwasP = gene_pval[[gene.name]]$gwasP
      }else{
        gene_POET_pval[[gene.name]]$snps = c(gene_pval[[gene.name]]$snps, gene.name)
        snps = gene_POET_pval[[gene.name]]$snps
        gene_POET_pval[[gene.name]]$zhat = c(gene_pval[[gene.name]]$zhat, gene_pval[[gene.name]]$rare_T)
        gene_POET_pval[[gene.name]]$gwasP = c(gene_pval[[gene.name]]$gwasP, gene_pval[[gene.name]]$rare_pval)
      }

      snps = gene_POET_pval[[gene.name]]$snps
      gene_POET_pval[[gene.name]]$LD <- corr.POET[snps, snps, drop = FALSE]
      LD_inverse = compute_ld_inverse(gene_POET_pval[[gene.name]][["LD"]], tol = 1e-5)
      gene_POET_pval[[gene.name]]$inv_LD = LD_inverse$inv_mat
      gene_POET_pval[[gene.name]]$eigenvalue = LD_inverse$ev
      gene_POET_pval[[gene.name]]$df = LD_inverse$df
      z_mat = as.matrix(gene_POET_pval[[gene.name]][["zhat"]])
      gene_POET_pval[[gene.name]]$chisq = as.numeric(t(z_mat) %*% gene_POET_pval[[gene.name]][["inv_LD"]] %*% z_mat)
      gene_POET_pval[[gene.name]]$pval = exp(pchisq(gene_POET_pval[[gene.name]][["chisq"]], df = gene_POET_pval[[gene.name]][["df"]], log.p = TRUE, lower.tail = FALSE))

      # gene-gene correlation
      gene1 = gene.name
      gene1.snps <- gene_POET_pval[[gene1]][["snps"]]

      st <- gene.loc$V3[as.character(gene.loc$V1)==gene1]
      ed <- gene.loc$V4[as.character(gene.loc$V1)==gene1]

      temp <- gene.loc.temp[-seq_len(i), , drop = FALSE]
      close.gene <- as.character(temp$V1[temp$V2 == chr & temp$V3 < ed + max.dist.bp & temp$V4 > st - max.dist.bp]) # 5MB
      close.gene <- close.gene[close.gene %in% sw.gene]
      close.gene <- close.gene[close.gene %in% names(gene_POET_pval)]

      if(length(close.gene)>0){
        gene_corr[[gene1]] = vector("list", length = length(close.gene))
        names(gene_corr[[gene1]]) = close.gene

        if(length(close.gene)>0){
          for(gene2 in close.gene){
            if(is.null(gene_pval[[gene2]]$rare_snps)){
              gene2.snps <- gene_pval[[gene2]][["snps"]]
            } else{
              gene2.snps <- c(gene_pval[[gene2]][["snps"]], gene2)
            }
            R <- corr.POET[c(gene1.snps, gene2.snps), c(gene1.snps, gene2.snps), drop = FALSE]

            p <- length(gene1.snps)
            q <- length(gene2.snps)
            gene1.invR.sqrt = compute_ld_inverse_sqrt(R[seq_len(p), seq_len(p), drop = FALSE], tol = 1e-5)
            gene2.invR.sqrt = compute_ld_inverse_sqrt(R[(p+1):(p+q), (p+1):(p+q), drop = FALSE], tol = 1e-5)

            three.parts = gene1.invR.sqrt %*% R[seq_len(p), (p+1):(p+q)] %*% gene2.invR.sqrt
            R.star <- diag(p+q)
            R.star[seq_len(p), (p+1):ncol(R.star)] <- three.parts
            R.star[(p+1):ncol(R.star), seq_len(p)] <- t(three.parts)

            R.star.chol <- cholesky_decomposition(R.star)
            cov.chisq <- compute_covariance_cholesky(Lmat = R.star.chol, p = p, q = q)
            corr.chisq <- cov.chisq / (sqrt(2*p) * sqrt(2*q))

            gene_corr[[gene1]][[gene2]]$R <- R
            gene_corr[[gene1]][[gene2]]$cov.chisq <- cov.chisq
            gene_corr[[gene1]][[gene2]]$corr.chisq <- corr.chisq

            # message("\t corr with ", gene2, "=", corr.chisq)
          }
        }
      }
    }
  }

  gene_corr.mat = diag(length(gene_corr))
  rownames(gene_corr.mat) <- names(gene_pval)
  colnames(gene_corr.mat) <- names(gene_pval)
  for (i in seq_len(length(gene_corr))) {
    if(i%%1000==1){
      message(i, "\t")
    }
    gene1 = names(gene_corr)[i]
    if(!is.null(gene_corr[[gene1]])){
      for (gene2 in names(gene_corr[[gene1]])) {
        gene_corr.mat[gene1, gene2] = gene_corr[[gene1]][[gene2]]$corr.chisq
        gene_corr.mat[gene2, gene1] = gene_corr[[gene1]][[gene2]]$corr.chisq
      }
    }
  }

  return(list(gene_POET_pval = gene_POET_pval, gene_corr.mat = gene_corr.mat))
}


compute_ld_inverse = function(mat, tol = 1e-8){
  if(!is.matrix(mat)){stop("Input has to be a matrix")}
  e = eigen(mat)
  U = e$vectors
  ev = e$values
  temp = try(ev > tol, silent=TRUE)
  if(is.character(temp)){
    ev <- Re(ev)
    U <- Re(U)
  }
  S = 1/ev*(ev>tol)
  S[is.na(S)] = 0
  inv_mat = U%*%diag(S)%*%t(U)
  df = sum(S>0)
  return(list(inv_mat = inv_mat, ev = ev, df = df))
}

compute_ld_inverse_sqrt = function(mat, tol = 1e-8){
  if(!is.matrix(mat)){stop("Input has to be a matrix")}
  e = eigen(mat)
  U = e$vectors
  ev = e$values
  temp = try(ev > tol, silent=TRUE)
  if(is.character(temp)){
    ev <- Re(ev)
    U <- Re(U)
  }
  S = 1/ev*(ev>tol)
  S[is.na(S)] = 0
  S.sqrt = sqrt(S)
  inv_sqrt_mat = U%*%diag(S.sqrt)%*%t(U)
  return(inv_sqrt_mat)
}

cholesky_decomposition <- function(R){
  R.chol = base::chol(R, pivot = TRUE)
  r <- attr(R.chol, 'rank')
  if(r < nrow(R)){
    R.chol[(r+1):nrow(R), (r+1):nrow(R)] <- 0
  }
  oo <- order(attr(R.chol, 'pivot'))
  R.pivot <- R.chol[,oo]

  return(t(R.pivot))
}

compute_covariance_cholesky <- function(Lmat, p, q){
  if(nrow(Lmat)!=(p+q) | ncol(Lmat)!=(p+q)){
    stop("Incorrect dimensions")
  }
  cov.chisq = 2 * sum(Lmat[(p+1):(p+q), seq_len(p)]^2)
  return(cov.chisq)
}


#' @title Burden test for rare variants
#'
#' @description We recover burden test statistics from GWAS summary statistics
#'  for the gene-level association testing of rare variants.
#'
#' @usage
#'  get_burden(genotype, gene_pval)
#'
#' @param genotype a genotype matrix
#' @param gene_pval a list of GWAS summary statistics for genes
#'
#' @return A list of gene-level testing for rare variants
#'
#' @examples
#' \dontrun{
#' rare_gene_pval.Demo <- get_burden(genotype = X_ref.Demo, gene_pval = rare_gene_pval.Demo)
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats
#' @export
get_burden <- function(genotype, gene_pval) {
  keep.snps <- unique(unlist(sapply(gene_pval, function(z){z$snps}), use.names = FALSE))
  keep.ref.idx <- which(rownames(genotype) %in% keep.snps)
  genotype <- genotype[keep.ref.idx, ]

  snpslength <- sapply(gene_pval, function(z){length(z$snps)})
  gene_pval <- gene_pval[which(snpslength <= 1000)]

  for (i in seq_len(length(gene_pval))) {
    gene.name = names(gene_pval)[i]
    if(i%%200==1){
      message("Burden test for ", i, ": ", gene.name)
    }

    snps <- gene_pval[[gene.name]][["snps"]]
    LD = make_ld_matrix_from_genotype(genotype, snps)
    rare.burden = burden_score(z = gene_pval[[gene.name]][["zhat"]], se = gene_pval[[gene.name]][["sehat"]], R = LD)
    gene_pval[[gene.name]]$U = rare.burden$U
    gene_pval[[gene.name]]$T = rare.burden$T
    gene_pval[[gene.name]]$bs = rare.burden$bs
    gene_pval[[gene.name]]$pval = exp(pchisq((gene_pval[[gene.name]][["T"]])^2, df = 1, log.p = TRUE, lower.tail = FALSE))
  }
  return(gene_pval)
}



make_ld_matrix_from_genotype <- function(genotype, all_snps) {
  # genotype: snps x individuals
  mat_dim <- length(all_snps)
  ld_matrix <- diag(mat_dim)

  # handles singleton genes
  if (mat_dim == 1){
    return(as.matrix(1))
  } else{
    idx = match(all_snps, rownames(genotype))
    X = t(genotype[idx, ])

    n <- nrow(X)
    Xs <- scale(X, center = TRUE, scale = TRUE)
    ld_matrix <- t(Xs) %*% Xs / (n-1)
    ld_matrix <- cov2cor(ld_matrix)
  }
  return(ld_matrix)
}

burden_score = function(z, se, R){
  U = z / se
  if(length(se)==1){
    diag_w = as.matrix(se^(-1))
  }else {
    diag_w = diag(se^(-1))
  }
  V = diag_w %*% R %*% diag_w
  nsnp = length(z)
  ksi <- as.vector(rep(1, nsnp))
  bs = as.numeric(t(ksi) %*% U)
  T = bs / sqrt(as.numeric(t(ksi) %*% V %*% ksi))
  return(list(U = U, V = V, T = T, bs = bs))
}




#' @title Combine gene-level testing for common and rare variants
#'
#' @description We combine gene-level chi-square association testing for common variants
#'  and burden test for rare variants. 1
#'
#' @usage
#'  get_combined(prunein_gene_pval, rare_gene_pval)
#'
#' @param prunein_gene_pval a list of gene-level testing for prune-in common variants
#' @param rare_gene_pval a list of gene-level testing for rare variants
#'
#' @return A list of combined gene-level testing
#'
#' @examples
#' \dontrun{
#' combine_gene_pval.Demo <- get_combined(prunein_gene_pval = prunein_gene_pval.Demo,
#'                                        rare_gene_pval = rare_gene_pval.Demo)
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @export
get_combined <- function(prunein_gene_pval, rare_gene_pval) {
  combine_gene_pval <- prunein_gene_pval

  for (i in seq_len(length(combine_gene_pval))) {
    gene.name = names(combine_gene_pval)[i]

    if(!is.null(rare_gene_pval[[gene.name]]) & length(rare_gene_pval[[gene.name]]$snps) <= 1000){
      combine_gene_pval[[gene.name]]$rare_snps = rare_gene_pval[[gene.name]]$snps
      combine_gene_pval[[gene.name]]$rare_betahat = rare_gene_pval[[gene.name]]$betahat
      combine_gene_pval[[gene.name]]$rare_sehat = rare_gene_pval[[gene.name]]$sehat
      combine_gene_pval[[gene.name]]$rare_zhat = rare_gene_pval[[gene.name]]$zhat
      combine_gene_pval[[gene.name]]$rare_gwasP = rare_gene_pval[[gene.name]]$gwasP
      combine_gene_pval[[gene.name]]$rare_U = rare_gene_pval[[gene.name]]$U
      combine_gene_pval[[gene.name]]$rare_T = rare_gene_pval[[gene.name]]$T
      combine_gene_pval[[gene.name]]$rare_bs = rare_gene_pval[[gene.name]]$bs
      combine_gene_pval[[gene.name]]$rare_pval = rare_gene_pval[[gene.name]]$pval

    }
  }
  return(combine_gene_pval)
}



#' @title Pseudo-SNP construction
#'
#' @description We collapse genotypes of all rare variants within a gene
#'  to construct a pseudo-SNP
#'
#' @usage
#'  construct_X_super(genotype, rare_gene_pval)
#'
#' @param genotype a genotype matrix
#' @param rare_gene_pval a list of GWAS summary statistics for rare variants
#'
#' @return A matrix of aggregated pseudo-SNP
#'
#' @examples
#' \dontrun{
#' X_super.Demo <- construct_X_super(genotype = X_ref.Demo, rare_gene_pval = rare_gene_pval.Demo)
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @export
construct_X_super <- function(genotype, rare_gene_pval) {
  snpslength <- sapply(rare_gene_pval, function(z){length(z$snps)})
  rare_gene_pval <- rare_gene_pval[which(snpslength <= 1000)]
  rare.snps <- unique(unlist(sapply(rare_gene_pval, function(z){z$snps}), use.names = FALSE))
  rare.ref.idx <- which(rownames(genotype) %in% rare.snps)
  genotype <- genotype[rare.ref.idx, ]

  X_super.list <- vector('list', length = length(rare_gene_pval))
  names(X_super.list) = names(rare_gene_pval)
  for (i in seq_len(length(rare_gene_pval))) {
    gene.name = names(rare_gene_pval)[i]
    if(i%%200==1){
      message("Construct pseudo-SNP for ", i, ": ", gene.name)
    }
    X_super.list[[gene.name]] <- construct_super_snps(genotype, rare_gene_pval[[gene.name]]$snps, super.name = gene.name)
  }
  X_super <- do.call(rbind, X_super.list)

  return(X_super)
}


construct_super_snps <- function(genotype, all_snps, super.name = NULL) {
  X <- genotype[all_snps, ,drop = FALSE]
  Gsup <- t(as.matrix(apply(X, 2, sum)))
  if(!is.null(super.name)){
    rownames(Gsup) <- super.name
  }
  return(Gsup)
}


#' @title Prioritizing trait-relevant tissue(s) and cell type(s)
#'
#' @description We devise a regression framework based on generalized least squares to identify risk loci enrichment
#'  to detect tissue- or cell-type-specific enrichment for a specific trait of interest.
#'
#' @usage
#'  prioritize_relevance(gene_pval, gene_corr.mat, gene_expr, gene.loc,
#'                       chrs = 1:22, verbose = FALSE)
#'
#' @param gene_pval a list of genes with gene-level association testing
#' @param gene_corr.mat a matrix of gene-gene correlations
#' @param gene_expr a data frame of gene expressions by tissues or cell types. For scRNA-seq data,
#'  this can be returned from \code{process_scRNA}.
#' @param gene.loc a data frame of gene dictionary, only including genes of interest in each transcriptomic profile
#'  (e.g., top 8000 highly variable genes from scRNA-seq or genes with great specificity score from GTEx v8).
#'  This can be returned from \code{process_scRNA}.
#' @param chrs chromosomes. Default is \code{1:22}.
#' @param verbose a logical value of diagnostic messages. Default is \code{FALSE}.
#'
#' @return EPIC enrichment results
#'
#' @examples
#' \dontrun{
#' gtex.enrichment.joint.Demo <- prioritize_relevance(gene_pval = combine_POET_gene_pval.Demo,
#'                                                    gene_corr.mat = combine_gene_corr.mat.Demo,
#'                                                    gene_expr = gtex.Demo, gene.loc = gtex.loc)
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats
#' @importFrom Matrix bdiag
#' @export
prioritize_relevance <- function(gene_pval, gene_corr.mat, gene_expr, gene.loc, chrs = 1:22, verbose = FALSE) {
  if(!all(names(gene_pval)==colnames(gene_corr.mat))){
    stop("Invalid dimension of gene-gene correlation!")
  }
  gene_pval <- gene_pval[names(gene_pval) %in% gene.loc$V1]
  gene_corr.mat <- gene_corr.mat[colnames(gene_corr.mat) %in% gene.loc$V1, colnames(gene_corr.mat) %in% gene.loc$V1]

  gene.loc <- gene.loc[order(gene.loc$V2, gene.loc$V3, gene.loc$V4), ]
  gene.keep = names(gene_pval)
  gene_expr.keep <- gene_expr[match(gene.keep, rownames(gene_expr)),]
  gene_expr.keep <- gene_expr.keep[match(gene.loc$V1, rownames(gene_expr.keep))[which(!is.na(match(gene.loc$V1, rownames(gene_expr.keep))))],]

  if(length(gene_pval) <= 1  | nrow(gene_corr.mat) <= 1 | nrow(gene_expr.keep) <= 1){
    stop("Invalid input of genes")
  }

  chisq <- sapply(gene_pval, function(z){z$chisq}, USE.NAMES = FALSE)
  df <- sapply(gene_pval, function(z){z$df}, USE.NAMES = FALSE)
  pval <- sapply(gene_pval, function(z){z$pval}, USE.NAMES = FALSE)
  pval.zero = which(pval<1e-300)
  chisq[pval.zero] <- qchisq(1e-300, df = df[pval.zero], lower.tail = FALSE)
  y <- chisq/df

  # remove the "Average" column
  tissues = colnames(gene_expr.keep)[-which(colnames(gene_expr.keep)=="Average")]

  chrs <- sort(unique(chrs))

  # Construct a block-diagonal matrix
  corr.list <- vector('list', length = length(chrs))
  chol.inv.list <- vector('list', length = length(chrs))
  y.list <- vector('list', length = length(chrs))
  y.reg.list <- vector('list', length = length(chrs))

  for (chr in chrs) {
    chr.idx <- which(chrs == chr)
    chr.genes <- gene.loc$V1[gene.loc$V2 == chr]
    block.idx <- match(chr.genes, gene.keep)[which(!is.na(match(chr.genes, gene.keep)))]
    corr.list[[chr.idx]] <- gene_corr.mat[block.idx, block.idx]

    y.list[[chr.idx]] <- y[block.idx]
    chol.inv.list[[chr.idx]] = compute_ld_inverse_sqrt(corr.list[[chr.idx]], tol = 1e-3)
    y.reg.list[[chr.idx]] <- chol.inv.list[[chr.idx]] %*% y.list[[chr.idx]]
  }
  y.reg <- unlist(y.reg.list[seq_len(length(chrs))])
  gene_expr.keep.reg <- Matrix::bdiag(chol.inv.list[seq_len(length(chrs))]) %*% as.matrix(gene_expr.keep)
  gene_expr.keep.reg <- as.matrix(gene_expr.keep.reg)

  enrichment_results <- rep(NA, length(tissues))
  names(enrichment_results) <- tissues

  for (i in tissues) {
    j = 1
    genes <- names(gene_pval)
    y.qc <- y.reg
    x1.qc <- gene_expr.keep.reg[,i]
    x2.qc <- gene_expr.keep.reg[,"Average"]
    df.qc <- df
    message("Prioritizing trait-relevant tissue(s) or cell type(s) for ", i)

    while(j <= 5){
      lm0 = lm(y.qc ~ x1.qc + x2.qc, weights = df.qc/2)
      cooksd <- cooks.distance(lm0)
      names(cooksd) <- genes
      cooksd_thres <- 0.2
      if(sum(cooksd > cooksd_thres) > 0){
        ig5 <- names(sort(cooksd, decreasing = TRUE)[1])
      }else{
        ig5 <- NULL
      }
      remove.idx <- (names(cooksd) %in% ig5)
      keep.idx <- !(names(cooksd) %in% ig5)

      if(verbose){
        message("\t Remove ", names(cooksd)[remove.idx])
      }
      lm <- lm(y.qc[keep.idx] ~ x1.qc[keep.idx] + x2.qc[keep.idx], weights = df.qc[keep.idx]/2)
      enrichment_results[i] <- exp(pt(summary(lm)$coefficients[2,3], df = summary(lm)$df[2], lower.tail = FALSE, log.p = TRUE))
      j = j + 1

      if(is.null(ig5)){
        break
      }

      genes <- genes[keep.idx]
      y.qc <- y.qc[keep.idx]
      x1.qc <- x1.qc[keep.idx]
      x2.qc <- x2.qc[keep.idx]
      df.qc <- df.qc[keep.idx]
    }
  }
  return(enrichment_results)
}







#' @title Identify a set of tissue- or cell-type-specific influential genes
#'
#' @description For a significantly enriched tissue or cell type,
#'  we further carry out a statistical influence test to identify a set of
#'  tissue- or cell-type-specific influential genes, using the DFBETAS statistics.
#' @usage
#'  influential_testing(gene_pval, gene_corr.mat, gene_expr, gene.loc,
#'                      chrs = 1:22, ct, verbose = FALSE)
#'
#' @param gene_pval a list of genes with gene-level association testing
#' @param gene_corr.mat a matrix of gene-gene correlations
#' @param gene_expr a data frame of gene expressions by tissues or cell types. For scRNA-seq data,
#'  this can be returned from \code{process_scRNA}.
#' @param gene.loc a data frame of gene dictionary, only including genes of interest in each transcriptomic profile
#'  (e.g., top 8000 highly variable genes from scRNA-seq or genes with great specificity score from GTEx v8).
#'  This can be returned from \code{process_scRNA}.
#' @param chrs chromosomes. Default is \code{1:22}.
#' @param ct a character of prioritized tissue or cell type
#' @param verbose a logical value of diagnostic messages. Default is \code{FALSE}.
#'
#' @return A list with components
#'  \item{influential.genes}{A list of tissue- or cell-type-specific influential genes}
#'  \item{p.dfbeta}{A ggplot object for DFBETA}
#'
#' @examples
#' \dontrun{
#' influential_test.Demo.liver <- influential_testing(gene_pval = combine_POET_gene_pval.Demo,
#'                                                    gene_corr.mat = combine_gene_corr.mat.Demo,
#'                                                    gene_expr = gtex.Demo, gene.loc = gtex.loc,
#'                                                    ct = "Liver")
#' influential.genes.liver = influential_test.Demo.liver$influential.genes
#' p.dfbeta.liver = influential_test.Demo.liver$p.dfbeta
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats olsrr ggplot2
#' @importFrom Matrix bdiag
#' @export
influential_testing <- function(gene_pval, gene_corr.mat, gene_expr, gene.loc, chrs = 1:22, ct, verbose = FALSE) {
  if(!all(names(gene_pval)==colnames(gene_corr.mat))){
    stop("Invalid dimension of gene-gene correlation!")
  }
  gene.loc <- gene.loc[order(gene.loc$V2, gene.loc$V3, gene.loc$V4), ]
  gene.keep = names(gene_pval)
  gene_expr.keep <- gene_expr[match(gene.keep, rownames(gene_expr)),]

  chisq <- sapply(gene_pval, function(z){z$chisq}, USE.NAMES = FALSE)
  df <- sapply(gene_pval, function(z){z$df}, USE.NAMES = FALSE)
  pval <- sapply(gene_pval, function(z){z$pval}, USE.NAMES = FALSE)
  pval.zero = which(pval<1e-300)
  chisq[pval.zero] <- qchisq(1e-300, df = df[pval.zero], lower.tail = FALSE)
  y <- chisq/df

  # remove the "Average" column
  tissues = colnames(gene_expr.keep)[-which(colnames(gene_expr.keep)=="Average")]

  chrs <- sort(unique(chrs))

  # Construct a block-diagonal matrix
  corr.list <- vector('list', length = length(chrs))
  chol.inv.list <- vector('list', length = length(chrs))
  y.list <- vector('list', length = length(chrs))
  y.reg.list <- vector('list', length = length(chrs))

  for (chr in chrs) {
    chr.idx <- which(chrs == chr)
    chr.genes <- gene.loc$V1[gene.loc$V2 == chr]
    block.idx <- match(chr.genes, gene.keep)[which(!is.na(match(chr.genes, gene.keep)))]
    corr.list[[chr.idx]] <- gene_corr.mat[block.idx, block.idx]

    y.list[[chr.idx]] <- y[block.idx]
    chol.inv.list[[chr.idx]] = compute_ld_inverse_sqrt(corr.list[[chr.idx]], tol = 1e-3)
    y.reg.list[[chr.idx]] <- chol.inv.list[[chr.idx]] %*% y.list[[chr.idx]]
  }
  y.reg <- unlist(y.reg.list[seq_len(length(chrs))])
  gene_expr.keep.reg <- Matrix::bdiag(chol.inv.list[seq_len(length(chrs))]) %*% as.matrix(gene_expr.keep)
  gene_expr.keep.reg <- as.matrix(gene_expr.keep.reg)


  j = 1
  genes <- names(gene_pval)
  y.qc <- y.reg
  x1.qc <- gene_expr.keep.reg[,ct]
  x2.qc <- gene_expr.keep.reg[,"Average"]
  df.qc <- df
  message("Prioritizing trait-relevant tissue(s) or cell type(s) for ", ct)

  while(j <= 5){
    lm0 = lm(y.qc ~ x1.qc + x2.qc, weights = df.qc/2)
    cooksd <- cooks.distance(lm0)
    names(cooksd) <- genes
    cooksd_thres <- 0.2
    if(sum(cooksd > cooksd_thres) > 0){
      ig5 <- names(sort(cooksd, decreasing = TRUE)[1])
    }else{
      ig5 <- NULL
    }
    remove.idx <- (names(cooksd) %in% ig5)
    keep.idx <- !(names(cooksd) %in% ig5)

    if(verbose){
      message("\t Remove ", names(cooksd)[remove.idx])
    }
    lm <- lm(y.qc[keep.idx] ~ x1.qc[keep.idx] + x2.qc[keep.idx], weights = df.qc[keep.idx]/2)
    j = j + 1

    if(is.null(ig5)){
      break
    }

    genes <- genes[keep.idx]
    y.qc <- y.qc[keep.idx]
    x1.qc <- x1.qc[keep.idx]
    x2.qc <- x2.qc[keep.idx]
    df.qc <- df.qc[keep.idx]
  }

  lm0 = lm(y.qc ~ x1.qc + x2.qc, weights = df.qc/2)

  # DFBETA
  obs <- NULL
  txt <- NULL
  dfb <- dfbetas(lm0)
  n <- nrow(dfb)
  np <- ncol(dfb)
  threshold <- 2/sqrt(n)
  myplots <- list()
  outliers <- list()
  dbetas <- dfb[, "x1.qc"]
  df_data <- data.frame(obs = seq_len(n), dbetas = dbetas)
  d <- ols_prep_dfbeta_data(df_data, threshold)
  d$txt <- ifelse(abs(d$dbetas)>threshold, genes, "")
  f <- ols_prep_dfbeta_outliers(d)
  p.dfbeta <- ggplot(d, aes(x = obs, y = dbetas, label = txt, ymin = 0, ymax = dbetas)) +
    geom_linerange() +
    geom_hline(yintercept = c(0, threshold, -threshold), colour = "red") +
    geom_point() + xlab("Genes") + ylab("DFBETAS") +
    ggtitle(paste("DFBETA for", ct)) +
    theme_bw() +
    theme(axis.text.x=element_blank(), text = element_text(size = 15)) +
    geom_text(hjust = -0.2, nudge_x = 0.15, size = 4,
              family = "serif", fontface = "italic",
              na.rm = TRUE) +
    annotate("text", x = Inf, y = Inf, hjust = 1.5, vjust = 2, size = 6,
             family = "serif", fontface = "italic", colour = "darkred",
             label = paste("Threshold:", round(threshold, 2)))
  influential.genes = d$txt[d$txt!=""]

  return(list(influential.genes = influential.genes, p.dfbeta = p.dfbeta))
}




if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("tissues", "pval"))
}

#' @title Plot tissue or cell type enrichment results
#'
#' @description Show barplots of inferred enrichment results for tissues or cell types
#'
#' @usage
#'  plot_relevance(enrichment_results)
#'
#' @param enrichment_results EPIC enrichment results, returned by \code{prioritize_relevance()}
#'
#' @return A ggplot object for enrichment visualization
#'
#' @examples
#' \dontrun{
#' enrichment_plot.Demo <- plot_relevance(gtex.enrichment.joint.Demo)
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import ggplot2
#' @export
plot_relevance <- function(enrichment_results) {
  o1 <- order(enrichment_results)
  enrichment_results = data.frame(tissues = names(enrichment_results), pval = unname(enrichment_results))

  p1 = ggplot(enrichment_results, aes(x = factor(tissues, levels = tissues[o1]), y = -log10(pval))) +
    geom_bar(stat="identity", position = position_dodge(), fill = "#92C5DE") +
    theme_bw() + xlab('') + ylab('-log10(p value)') +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 55, hjust = 1, size = 10)) +
    geom_hline(yintercept = -log10(0.05/nrow(enrichment_results)), linetype = 2)

  return(p1)
}








