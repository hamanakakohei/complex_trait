#!/usr/bin/env bash
#
# 目的：
# molQTLとGWASのcolocalizeを示す用
# なのでeGeneとかを指定する仕組みになっている
# gwas sumstaを読むときのフォーマットは以下のファイルに自作のを追記する：
# --> ~/miniconda3/envs/gwaslab3/lib/python3.10/site-packages/gwaslab/data/formatbook.json
#
# 入力
# - locuscompare_input.txt
# -- タブ区切りで以下の列順
# -- rsid: 図示したいバリアント（eQTLリードとか）
# -- chr: 上のchr、これでフィルターしないとsumsta読むのがきつい
# -- phenotype_id: QTLのeGeneとか、これでtensorqtl nominal結果をフィルターして図示する
# -- f: gwas sumstaファイル名、ディレクトリはslurmスクリプト内で勝手に指定している、、、
# -- fmt: 上のをgwaslabで読むときのフォーマット (自作したものでもOK)
# -- popu: locuscompareで図示するときのLD pop (UKB, EAS, EUR, AFR, etc.)
# -- 以降の列は無関係
# - rnaseq430.nominal.cis_qtl_pairs.chr
# -- tensorqtl nominalの結果のparquetファイルのchr番号前までのprefix
#
# 出力
# locuscompare画像
set -euo pipefail


# 1
STA=34
END=$(wc -l < inputs/uniq_sumsta.list)
END=34

sbatch \
  --array=${STA}-${END}%100 \
  slurm/pipeline_qtl_gwas_coloc_array.01.slurm \
  inputs/uniq_sumsta.list \
  inputs/tensorqtl_nominal/rnaseq430.nominal.cis_qtl_pairs.chr
  #--array=2,11,61,70,71,75,85,145,146,147,163,164,179%50 \
  #--array=4,69,80,127,128%50 \
  #--array=9,15,19,62-66,72-74,100-103,108-109,113,124,148-159,161,180,181%50 \
  #inputs/locuscompare_input.txt \
  #--cpus-per-task=5 \


  #58,59,91,262,277 --> vcfなし
  #142 --> KeyError
  #75,76,77,etc --> Extreme INFO/RS value --> see https://github.com/Cloufield/gwaslab/issues/191
  #101,102,103,etc --> pd.Categorical
  #18,86,96,etc(geneatlas系) --> Index(['rsID', 'CHR', etc.で終わるがエラーではない
  #
  #
  #
  #
  #
  # 34315874-GCST90013664-EFO_0004736.h.tsv.gz: BETA, OR, SEなしだが、ＮがわかっているからOK???
  #
  # 以下、quantなのにMAFない or EAFかどうか不明なのでsdYを1とみなしてcolocする
  # 34315874-GCST90013664-EFO_0004736.h.tsv.gz
  # GCST90019494.h.tsv.gz
  # GCST90019498.h.tsv.gz
  # GCST90019502.h.tsv.gz
  # GCST90019504.h.tsv.gz
  # GCST90019506.h.tsv.gz
  # GCST90019509.h.tsv.gz
  # GCST90019515.h.tsv.gz
  # GCST90019518.h.tsv.gz
  # GCST90019522.h.tsv.gz
  # GCST90019525.h.tsv.gz
  # imputed.allWhites.1697-0.0.chr
  # imputed.allWhites.1727-0.0.chr
  # imputed.allWhites.20015-0.0.chr
  # imputed.allWhites.23106-0.0.chr
  # imputed.allWhites.30010-0.0.chr
  # imputed.allWhites.30120-0.0.chr
  # imputed.allWhites.30140-0.0.chr
  # imputed.allWhites.50-0.0.chr
