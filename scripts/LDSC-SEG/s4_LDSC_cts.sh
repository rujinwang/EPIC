#!/bin/bash

trait=blah_trait_blah

ldsc_dir=/nas/longleaf/home/rujin/ldsc
pine=/pine/scr/r/u/rujin
out_dir=${pine}/EPIC/gtex_v8/ldsc_results/cts/${trait}
bim_dir=${pine}/EPIC/tutorial/ldsc/1000G_EUR_Phase3_plink

[ ! -d $out_dir ] && mkdir -p $out_dir

#-------------
# calculate tissue specific
#-------------

# create a file GTEx_v8.ldcts!!!!!!

cd ${pine}/EPIC/gtex_v8/ldsc_results

echo "$trait"
python $ldsc_dir/ldsc.py \
	--h2-cts $pine/EPIC/ldsc_input/$trait/ldsc_${trait}.sumstats.gz \
	--ref-ld-chr ${pine}/EPIC/tutorial/ldsc/1000G_EUR_Phase3_baseline/baseline. \
	--out $out_dir/${trait}_GTEx_v8 \
	--ref-ld-chr-cts ${pine}/EPIC/gtex_v8/ldsc_results/GTEx_v8.ldcts \
	--w-ld-chr ${pine}/EPIC/tutorial/ldsc/weights_hm3_no_hla/weights.

