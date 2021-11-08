#!/bin/bash
trait=blah_trait_blah

dir=/pine/scr/r/u/rujin/EPIC

if [ $trait == T2Db ]; then
	input_dir=$dir/sc_pancreas/magma_input/$trait
else
	input_dir=$dir/sonawane/magma_input/$trait
fi

out_dir=$dir/gtex_v8/magma_v1.08_results/$trait
magma_dir=/nas/longleaf/home/rujin/magma_v1.08_static
ref_dir=/pine/scr/r/u/rujin/scRNAseqGWAS/MAGMA

if [ ! -e $input_dir/magma_${trait}_snps.txt ]; then
	awk '{print $3 "\t" $1 "\t" $2}' $input_dir/${trait}_magma_gwas_input.txt | tail -n +2 > $input_dir/magma_${trait}_snps.txt
fi

[ ! -d $out_dir ] && mkdir -p $out_dir

# annotate
$magma_dir/magma --annotate window=10,1.5 \
				--snp-loc "$input_dir/magma_${trait}_snps.txt" \
				--gene-loc "$dir/gtex_v8/gtex_v8.gene.noMHC.loc" \
				--out "$out_dir/${trait}"				
				
$magma_dir/magma --bfile "$ref_dir/g1000_eur/g1000_eur" \
	--pval "$input_dir/${trait}_magma_gwas_input.txt" use=rsid,P ncol=N \
	--gene-annot "$out_dir/${trait}.genes.annot" \
	--out "$out_dir/${trait}"
	
# Corrections for covariates

# GTEx v8	
$magma_dir/magma --gene-results "$out_dir/${trait}.genes.raw" \
	--gene-covar "$dir/gtex_v8/log2_tpm_median.txt" --model direction=pos condition='Average' \
	--out "$out_dir/${trait}_cond"




				