#!/bin/bash

#Test for the firth batch
plink2 --pfile /vol/projects/CIIM/Lyme_GWAS/GWAS/Genetics/output/batch1/PostImp/merge_filtered_rsid --pheno iid-only /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/covars/pheno.tsv --covar iid-only /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/covars/covariates.tsv --covar-name sex --glm --ci 0.95 --out /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_500FG

#Test for the second batch
plink2 --pfile /vol/projects/CIIM/Lyme_GWAS/GWAS/Genetics/output/batch2/PostImp/merge_filtered_rsid --pheno iid-only /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/covars/pheno.tsv --covar iid-only /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/covars/covariates.tsv --covar-name sex --glm --ci 0.95 --out /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme2_300BCG

#Metal for both
metal metal.sh
mv METAANALYSIS1.TBL.info /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_lyme2.metal.tbl.info
mv METAANALYSIS1.TBL /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_lyme2.metal.tbl

#Metal for both and finngen
#First, we need to cut only the rsid
cat  /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_500FG.pheno.glm.logistic.hybrid | sed 's/ID/%ID/g'  | cut -f 2 -d % > /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_500FG_rsid.pheno.glm.logistic.hybrid

cat  /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme2_300BCG.pheno.glm.logistic.hybrid | sed 's/ID/%ID/g'  | cut -f 2 -d % > /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme2_300BCG_rsid.pheno.glm.logistic.hybrid

metal metal_finngen.sh
mv METAANALYSIS1.TBL.info /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_lyme2_finngen.metal.tbl.info
mv METAANALYSIS1.TBL /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_lyme2_finngen.metal.tbl
