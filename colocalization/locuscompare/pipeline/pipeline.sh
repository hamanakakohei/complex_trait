#!/usr/bin/env bash
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate gwaslab3


LEAD_RS=""
CHR=""
GWA_SUMSTA=""
GWA_FMT=""
QTL_SUMSTA=""
ENSG=""
LD_POP=""
PLOT_WINDOW="1000"
THREAD=5
OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --lead_rs)     LEAD_RS="$2";    shift 2;;
    --chr)         CHR="$2";        shift 2;;
    --gwa_sumsta)  GWA_SUMSTA="$2"; shift 2;;
    --gwa_fmt)     GWA_FMT="$2";    shift 2;;
    --qtl_sumsta)  QTL_SUMSTA="$2"; shift 2;;
    --ensg)        ENSG="$2";       shift 2;;
    --plot_window) WINDOW="$2";     shift 2;;
    --ld_pop)      LD_POP="$2";     shift 2;;
    --thread)      THREAD="$2";     shift 2;;
    --out)         OUT="$2";        shift 2;;
    *) echo "Unknown option: $1";  exit 1;;
  esac
done

GWA_SUMSTA_SUBSET=results/01/$(basename $GWA_SUMSTA).$LEAD_RS.gz
QTL_SUMSTA_SUBSET=results/02/$(basename $QTL_SUMSTA).$ENSG.txt.gz
GWA_SUMSTA_ORIG=$GWA_SUMSTA


# 0: gwas sumstaをgwaslab用にmunge
# 方針
# - formatbookで対応するのが一番シンプルで柔軟
# -- ~/miniconda3/envs/gwaslab3/lib/python3.10/site-packages/gwaslab/data/formatbook.json
# - それでも無理なら事前にファイルを以下のようにいじる
# - 結構ファイルが重いので事前に染色体で削るのがバランスが良い
TEMP1=$(mktemp)
TEMP2=$(mktemp)   # pan-ukb 用
trap 'rm "$TEMP1" "$TEMP2"' EXIT

# ---- フォーマット別フィルタ ----
case "$GWA_FMT" in
  gwascatalog_hm)
    zcat "$GWA_SUMSTA" \
      | awk -v CHR="$CHR" 'NR==1 || $3==CHR' \
      > "$TEMP1"
    ;;

  fastgwa|ssf)
    zcat "$GWA_SUMSTA" \
      | awk -v CHR="$CHR" 'NR==1 || $1==CHR' \
      | sed -E 's/rs756813921([[:space:]]|$)/rs34498239\1/g' \
      > "$TEMP1"
    ;;

  vcf)
    # ヘッダーを飛ばして、列名を直す
    zcat "$GWA_SUMSTA" \
      | awk -v CHR="$CHR" '
          BEGIN{OFS="\t"}
          /^#CHROM/ { $10="Study_1"; print }
          !/^#/ && ($1==CHR || $1=="chr"CHR){ print }
        ' \
      > "$TEMP1"
    ;;

  geneatlas)
    # タブ区切り+列名を直す
    # SNP列はrsIDとSNPIDが混ざっていて特定の染色体を抜き出せない
    zcat ${GWA_SUMSTA}${CHR}.csv.gz \
      | awk '
          NR==1{print "SNP ALLELE iscores NBETA NSE PV"}
          NR>1 {print}
        ' \
      | tr ' ' '\t' \
      > "$TEMP1"
    ;;

  pan-ukb_*)
    # gwaslabが.bgzを読めない
    zcat "$GWA_SUMSTA" \
      | awk -v CHR="$CHR" 'NR==1 || $1==CHR' \
      > "$TEMP1"

    # OOMならないように列を削る
    BAD_COLS="AFR|AMR|CSA|EAS|MID|heterogeneity|low_confidence|af_c" # EUR
    COLS=$(head -n1 "$TEMP1" \
             | tr '\t' '\n' \
             | nl -w1 -s: \
             | grep -vE "$BAD_COLS" \
             | cut -d: -f1 \
             | paste -sd, -)

    cut -f$COLS $TEMP1 > $TEMP2
    TEMP1=$TEMP2
    ;;

  NealeLab)
    # gwaslabが.bgzを読めない
    zcat $GWA_SUMSTA > $TEMP1
    ;;
  *)
    # tempファイル名から圧縮の有無を推測できないので解凍しておく
    ( [[ $GWA_SUMSTA == *.gz ]] && zcat $GWA_SUMSTA || cat $GWA_SUMSTA ) > $TEMP1
    ;;
esac

# ---- 最終的な出力を GWA_SUMSTA として置き換える ----
GWA_SUMSTA=$TEMP1


## 1: gwas sumstaで必要な部分を取る
#python scripts/01_subset_gwas_sumsta_interval.py \
#  --sumsta $GWA_SUMSTA \
#  --fmt $GWA_FMT \
#  --rsid $LEAD_RS \
#  --window $WINDOW \
#  --threads $THREAD \
#  --out $GWA_SUMSTA_SUBSET \
#  #> logs/1.$ANALYSIS.log 2>&1
#
#
## 2: qtl sumstaで必要な部分を取る
#python scripts/02_subset_qtl_sumsta_by_phen.py \
#  --qtl_nominal $QTL_SUMSTA \
#  --phen $ENSG \
#  --rsid $LEAD_RS \
#  --window $WINDOW \
#  --threads $THREAD \
#  --out $QTL_SUMSTA_SUBSET \
#  #> logs/2.$ANALYSIS.log 2>&1
#
#
## 3: locuscompareRでプロットする
#conda activate locuscompare
#Rscript scripts/03_run_locuscompare.R \
#  --sumsta1 $GWA_SUMSTA_SUBSET \
#  --sumsta2 $QTL_SUMSTA_SUBSET \
#  --title1 gwas \
#  --title2 qtl \
#  --rsid $LEAD_RS \
#  --pop $LD_POP \
#  --assembly hg38 \
#  --out $OUT \
#  #> logs/3.$ANALYSIS.log 2>&1


# 4: gwaslabでgtfとプロットする
# to do: 中身が01とまろかぶりなので、その部分をsumstaの共通関数とする
GTF=inputs/xxx.gtf

python scripts/04_regional_plot_by_gwaslab.py \
  --sumsta $GWA_SUMSTA \
  --fmt $GWA_FMT \
  --rsid $LEAD_RS \
  --gtf $GTF \
  --pop $LD_POP \
  --window $WINDOW \
  --threads $THREAD \
  --out results/04/$LEAD_RS.$ENSG.$(basename $GWA_SUMSTA_ORIG).$LD_POP.gwaslab.png
