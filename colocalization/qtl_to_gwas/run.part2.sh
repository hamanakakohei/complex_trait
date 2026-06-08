#!/usr/bin/env bash
#
# 目的： eqtl - gwas coloc
#
# 入力
#
# 出力
set -euo pipefail


# 1. gwasスタディのさらなるメタ情報をcausaldb downstream analysisファイルにくっつける
eval "$(conda shell.bash hook)"
conda activate misc_20250301

python - <<'EOF'
import pandas as pd
pd.read_table("inputs/part2/downstream_analysis_input.txt").merge(
    pd.read_table("results/part1/02/gwas_info.txt"),
    on="f"
).to_csv("results/part2/01/downstream_analysis_input.AddGwasInfo.txt", sep="\t", index=False, na_rep="NA")
EOF


# 2.
INPUT=results/part2/01/downstream_analysis_input.AddGwasInfo.txt
ASSOC_PREFIX=inputs/part2/tensorqtl_nominal/rnaseq430.nominal.cis_qtl_pairs.
GENO_PREFIX=inputs/part2/plink_bed/ALL.correctRefAlt.norm.
PHEN=inputs/part2/rnaseq430.expression.bed.gz
COV=inputs/part2/rnaseq430.combined_covariates.txt

STA=2 # 1行目はヘッダー
END=$(wc -l < $INPUT)
#END=57

sbatch \
        --array=${STA}-${END}%100 \
        slurm/pipeline_qtl_gwas_coloc_array.02.slurm \
        $INPUT \
        $ASSOC_PREFIX \
        $GENO_PREFIX \
        $PHEN \
        $COV
        #--array=2-4,25,26,33,69,87,91,92,118,230-231,240-241,249,251-253,255%50 \
