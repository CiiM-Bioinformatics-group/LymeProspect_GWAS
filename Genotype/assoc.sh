#!/bin/bash

#Test for the firth batch
plink2 --pfile output/batch1/PostImp/merge_filtered_rsid --pheno iid-only covars/pheno.tsv --covar iid-only covars/covariates.tsv --covar-name sex --glm --ci 0.95 --out out/lyme1_500FG

#Test for the second batch
plink2 --pfile output/batch2/PostImp/merge_filtered_rsid --pheno iid-only /covars/pheno.tsv --covar iid-only covars/covariates.tsv --covar-name sex --glm --ci 0.95 --out out/lyme2_300BCG

#Metal for both
metal metal.sh
mv METAANALYSIS1.TBL.info out/lyme1_lyme2.metal.tbl.info
mv METAANALYSIS1.TBL out/lyme1_lyme2.metal.tbl

#Metal for both and finngen
#First, we need to cut only the rsid
cat  out/lyme1_500FG.pheno.glm.logistic.hybrid | sed 's/ID/%ID/g'  | cut -f 2 -d % > out/lyme1_500FG_rsid.pheno.glm.logistic.hybrid

cat  out/lyme2_300BCG.pheno.glm.logistic.hybrid | sed 's/ID/%ID/g'  | cut -f 2 -d % > out/lyme2_300BCG_rsid.pheno.glm.logistic.hybrid

metal metal_finngen.sh
mv METAANALYSIS1.TBL.info out/lyme1_lyme2_finngen.metal.tbl.info
mv METAANALYSIS1.TBL out/lyme1_lyme2_finngen.metal.tbl
