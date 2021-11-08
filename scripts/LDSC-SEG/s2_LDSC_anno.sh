#!/bin/bash


#-------------
# make annotation
#-------------

ldsc_dir=/nas/longleaf/home/rujin/ldsc
out_dir=/pine/scr/r/u/rujin/EPIC/gtex_v8/ldsc_results
bim_dir=/pine/scr/r/u/rujin/EPIC/tutorial/ldsc/1000G_EUR_Phase3_plink

[ ! -d $out_dir/anno ] && mkdir -p $out_dir/anno

tissue=blah_tissue_blah


if [ $tissue == control ]; then
	# Generate control gene set
	# The control gene sets are all genes that were included in the analysis to
	# determine differential expression. 

	part_start=1
	part_end=22
	for chr in `seq $part_start $part_end`; do
		python $ldsc_dir/make_annot.py \
			--gene-set-file $out_dir/SEG/control.txt \
			--gene-coord-file /pine/scr/r/u/rujin/scRNAseqGWAS/ldsc/example_ldsc/ENSG_coord.txt \
			--windowsize 100000 \
			--bimfile $bim_dir/1000G.EUR.QC.${chr}.bim \
			--annot-file $out_dir/anno/control_chr${chr}.annot.gz
	done
else
	part_start=1
	part_end=22
	for chr in `seq $part_start $part_end`; do
		python $ldsc_dir/make_annot.py \
			--gene-set-file $out_dir/SEG/top10pct_GTEx_v8_${tissue}.txt \
			--gene-coord-file /pine/scr/r/u/rujin/scRNAseqGWAS/ldsc/example_ldsc/ENSG_coord.txt \
			--windowsize 100000 \
			--bimfile $bim_dir/1000G.EUR.QC.${chr}.bim \
			--annot-file $out_dir/anno/top10pct_GTEx_v8_${tissue}_chr${chr}.annot.gz
	done
fi







