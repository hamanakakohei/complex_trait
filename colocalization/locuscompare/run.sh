#!/usr/bin/env bash
#
# 目的：
# molQTLとGWASのcolocalizeを示す用
# なのでeGeneとかを指定する仕組みになっている
# gwas sumstaを読むときのフォーマットは以下のファイルに自作のを追記する：
# --> ~/miniconda3/envs/gwaslab3/lib/python3.10/site-packages/gwaslab/data/formatbook.json
#
# 入力：
# 1. locuscompare_input.txt
# タブ区切りで以下の列順
# - rsid: 図示したいバリアント（eQTLリードとか）
# - chr: 上のchr、これでフィルターしないとsumsta読むのがきつい
# - phenotype_id: QTLのeGeneとか、これでtensorqtl nominal結果をフィルターして図示する
# - f: gwas sumstaファイル名、ディレクトリはslurmスクリプト内で勝手に指定している、、、
# - fmt: 上のをgwaslabで読むときのフォーマット (自作したものでもOK)
# - popu: locuscompareで図示するときのLD pop (UKB, EAS, EUR, AFR, etc.)
# - 以降の列は無関係
#
# 2. rnaseq430.nominal.cis_qtl_pairs.chr
# tensorqtl nominalの結果のparquetファイルのchr番号前までのprefix
#
# 出力
# locuscompare画像
set -euo pipefail

STA=1
END=$(wc -l < inputs/locuscompare_input.txt)
END=$(expr $END - 1)
#END=1

sbatch \
  --cpus-per-task=12 \
  --array=${STA}-${END}%20 \
  slurm/pipeline_array.slurm \
  inputs/locuscompare_input.txt \
  inputs/qtl_sumsta/rnaseq430.nominal.cis_qtl_pairs.chr
  #--array=1,17,85,95,96,103,104,108,111,121,182,183,184,236,237,238,255,256,273,288%20 \
