#!/usr/bin/env bash
#
# 0：あるeQTLがそのeGene以外の遺伝子との有意な関連（nominal P）とPIP値をまとめた表を出す
# - 入力
# -- 1. 着目しているeQTLのphenotype_id列、variant_id列（と他の列）からなるファイル
# -- 2. tensorqtl all signif pairs結果ファイル
# -- 3. tensorqtl susie結果ファイル
# - 出力
# -- 入力ファイル１に他の遺伝子とのnominal PとPIPを結合したもの
#
# １：複数の遺伝子のeqtl間でcolocする、ここでは0で有意な関連があった別遺伝子とeGene間
# - 入力
# -- 1. tensorqtl nominal結果（の染色体番号前までのprefix）
# -- 2. tensorqtl解析時に用いたgenotypeのplink bed/bim/famファイルのprefix
# -- 3. tensorqtl解析時に用いたPHENファイル
# -- 4. tensorqtl解析時に用いたCOVファイル
# - 出力
# -- 各遺伝子におけるsusie結果、遺伝子間のcoloc (susie or abf)結果、etc.
# -- ざっくりregional plot画像
#
# To do
# - 色んなオプションが01_qtl_qtl_coloc.slurm内で指定するようになってしまっている、このスクリプト内に出してくる
set -euo pipefail


# 0.
ALL_SIGNIF=inputs/all_signif_pairs.txt.gz
CS95=inputs/rnaseq430.SuSiE_summary.parquet

python scripts/00.py \
  --new_gene_qtl inputs/trait-associated_new_gene_qtl.txt \
  --tensorqtl_all_signif_res $ALL_SIGNIF \
  --tensorqtl_susie_res $CS95 \
  --out results/00/trait-associated_new_gene_qtl.summary.txt


# 1.
phen1_lead_rank_phen2=results/00/trait-associated_new_gene_qtl.summary.txt
ASSOC_PREFIX=inputs/tensorqtl_nominal/rnaseq430.nominal.cis_qtl_pairs.
GENO_PREFIX=inputs/plink_bed/ALL.correctRefAlt.norm.
PHEN=inputs/rnaseq430.expression.bed.gz
COV=inputs/rnaseq430.combined_covariates.txt

STA=2 # 1行目はヘッダー
END=$(wc -l < $phen1_lead_rank_phen2)

sbatch --array=${STA}-${END}%40 \
        slurm/01_qtl_qtl_coloc.slurm \
        $phen1_lead_rank_phen2 \
        $ASSOC_PREFIX \
        $GENO_PREFIX \
        $PHEN \
        $COV
