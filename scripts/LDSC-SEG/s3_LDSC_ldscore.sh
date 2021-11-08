#!/bin/bash



ldsc_dir=/nas/longleaf/home/rujin/ldsc
out_dir=/pine/scr/r/u/rujin/EPIC/gtex_v8/ldsc_results
bim_dir=/pine/scr/r/u/rujin/EPIC/tutorial/ldsc/1000G_EUR_Phase3_plink

[ ! -d $out_dir/ldscore ] && mkdir -p $out_dir/ldscore
#-------------
# calculate ldscore
#-------------

chr=blah_chr_blah
tissue=blah_tissue_blah
echo ${chr}
echo $tissue

if [ $tissue == control ]; then
	# control 
	python $ldsc_dir/ldsc.py \
		--l2 \
		--bfile $bim_dir/1000G.EUR.QC.${chr} \
		--ld-wind-cm 1 \
		--annot $out_dir/anno/control_chr${chr}.annot.gz \
		--thin-annot \
		--out $out_dir/ldscore/control_chr${chr}
else
	python $ldsc_dir/ldsc.py \
		--l2 \
		--bfile $bim_dir/1000G.EUR.QC.${chr} \
		--ld-wind-cm 1 \
		--annot $out_dir/anno/top10pct_GTEx_v8_${tissue}_chr${chr}.annot.gz \
		--thin-annot \
		--out $out_dir/ldscore/top10pct_GTEx_v8_${tissue}_chr${chr}
fi
