library('biomaRt')

pval_ACAT=function(pvals){
  # Cauchy combination
  Pvals=pvals
  Pvals[Pvals>0.999]=0.999
  is.small<-(Pvals<1e-15)
  Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
  Pvals[is.small]<-1/Pvals[is.small]/pi
  cct.stat<-mean(Pvals)
  pval.cauchy<-pcauchy(cct.stat,lower.tail = F)
  return(pval.cauchy)
}

for(trait in c('SCZ', 'BIP','SCZBIP','HDL','LDL','TC','TG','T2Db')){
  
  cat(trait,'\n\n\n')  
  # trait GWAS
  load(paste0('utils/',trait,'_aligned.rda'))
  head(gwas.align) # SNP-level testing statistics
  hist(gwas.align$P)
  
  # snp to gene assignment
  load(paste0('utils/',trait,'_combine_0.01_QC.rda'))
  combine_gene_pval[[1]]
  
  # EPIC combined
  load(paste0('combined/gtex/',trait,'_combine_0.01_chisq_iPOET_QC.noMHC.rda'))
  length(combine_POET_gene_pval)
  
  # EPIC common calling
  load(paste0('common/gtex_v8/',trait,'_common_0.01_chisq_POET_QC.noMHC.rda'))
  length(prunein_POET_gene_pval)
  
  if(!all(names(combine_POET_gene_pval)==names(prunein_POET_gene_pval))) break
  
  # EPIC rare calling
  load(paste0('rare/',trait,'_rare_0.01_burden.rda'))
  length(rare_gene_pval)
  
  
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id",
                  attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values=names(prunein_POET_gene_pval),mart= mart)
  G_list=as.matrix(G_list)
  
  prunein_POET_gene_pval=prunein_POET_gene_pval[match(G_list[,1], names(prunein_POET_gene_pval))]
  combine_POET_gene_pval=combine_POET_gene_pval[match(G_list[,1], names(combine_POET_gene_pval))]
  
  all(G_list[,1]==names(prunein_POET_gene_pval))
  all(G_list[,1]==names(combine_POET_gene_pval))
  
  gene.i=which(G_list[,2]=='PPARG') # T2Db
  gene.i=which(G_list[,2]=='APOB') # LDL
  gene.i=which(G_list[,2]=='APOE') # LDL
  
  
  nCommon=nRare=ACAT_joint=ACAT_common=ACAT_rare=EPIC_joint=EPIC_common=EPIC_rare=rep(NA,nrow(G_list))
  for(i in 1:nrow(G_list)){
    gene.i=G_list[i,1]
    
    nCommon[i]=length(prunein_POET_gene_pval[[gene.i]]$snps)

    
    # Common
    EPIC_common[i]=prunein_POET_gene_pval[[gene.i]]$pval # common variants p-value
    length(prunein_POET_gene_pval[[gene.i]]$gwasP) # number of common SNPs
    ACAT_common[i]=pval_ACAT(prunein_POET_gene_pval[[gene.i]]$gwasP)
    
    # Combined
    EPIC_joint[i]=combine_POET_gene_pval[[gene.i]]$pval
    
    if(!is.null(rare_gene_pval[[gene.i]])){
      # ACAT_combined
      ACAT_joint[i]=pval_ACAT(c(prunein_POET_gene_pval[[gene.i]]$gwasP, rare_gene_pval[[gene.i]]$gwasP))
      
      # Rare
      EPIC_rare[i]=rare_gene_pval[[gene.i]]$pval # rare variants p-value
      length(rare_gene_pval[[gene.i]]$gwasP) # number of rare SNPs
      ACAT_rare[i]=pval_ACAT(rare_gene_pval[[gene.i]]$gwasP)
      
      nRare[i]=length(rare_gene_pval[[gene.i]]$snps)
    } else{
      ACAT_joint[i]=ACAT_common[i]
    }
   
  }
  
  output=data.frame(G_list, nCommon, nRare, EPIC_joint, EPIC_common, EPIC_rare,
                    ACAT_joint, ACAT_common, ACAT_rare)
  
  
  # MAGMA calling results
  magma=read.table(paste0('magma/',trait,'.genes.out'), header = T)
  head(magma)
  
  
  # SKAT calling results
  load(paste0('SKAT/Gene_pvalue_',trait,'_0.01_full_results.rda'))
  head(full_gene_level_results_matrix)
  skat=full_gene_level_results_matrix[,c('gene','SKAT_joint','SKAT_common','SKAT_rare')]
  
  
  genes=intersect(intersect(output[,1], magma[,1]), skat[,1])
  
  skat=skat[match(genes, skat[,1]),-1]
  magma=magma[match(genes, magma[,1]),]
  output=output[match(genes, output[,1]),]
  output=cbind(output,skat,MAGMA=magma$P)
  
  # output[which(output$hgnc_symbol=='CNTN4'),]
  
  save(output, file=paste0(trait,'_common_rare.rda'))
  
}

# For T2Db
output[which(output$hgnc_symbol=='SCD5'),]
output[which(output$hgnc_symbol=='HNF4A'),]
