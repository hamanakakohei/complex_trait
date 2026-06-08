#!/usr/bin/env bash
#
# 注意点：
# - 50.v1.1.fastGWA.gzのrsIDが古くてrs756813921 -> rs34498239 に変えている
# - 本当はgwaslab内でrsIDを特定のdbSNPバージョンに揃えたらよいが、、、
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

GWA_SUMSTA_BASE=$(basename $GWA_SUMSTA)
GWA_SUMSTA_SUBSET=results/01/$GWA_SUMSTA_BASE.chr$CHR.$LEAD_RS.gz
#QTL_SUMSTA_SUBSET=results/02/$(basename $QTL_SUMSTA).$ENSG.txt.gz


# 0: gwas sumstaをgwaslab用にmunge
# 方針
# - formatbookで対応するのが一番シンプルで柔軟
# -- ~/miniconda3/envs/gwaslab3/lib/python3.10/site-packages/gwaslab/data/formatbook.json
# - それでも無理なら事前にファイルを以下のようにいじるがここでは最低限
# - 結構ファイルが重いので事前に染色体で削るのがバランスが良い
# - でgwaslabで各フォーマットに合わせてしっかりいじる
TEMP1=$(mktemp)
TEMP2=$(mktemp)
TEMP3=$(mktemp)
trap 'rm $TEMP1 $TEMP2 $TEMP3' EXIT

# ---- フォーマット別フィルタ ----
case "$GWA_FMT" in
  gwascatalog_hm)
    zcat "$GWA_SUMSTA" \
      | awk -v CHR="$CHR" 'NR==1 || $3==CHR' \
      > "$TEMP1"
    ;;

  ssf)
    zcat "$GWA_SUMSTA" \
      | awk -v CHR="$CHR" 'NR==1 || $1==CHR' \
      | sed -E 's/rs756813921([[:space:]]|$)/rs34498239\1/g' \
      > "$TEMP1"
    ;;

  #fastgwa_v1.0)
  #  zcat "$GWA_SUMSTA" \
  #    | awk -v CHR="$CHR" 'NR==1 || $1==CHR' \
  #    > "$TEMP1"
  #  ;;

  fastgwa)
    zcat "$GWA_SUMSTA" \
      | awk -v CHR="$CHR" -v OFS="\t" '
      NR==1 {
        for (i=1; i<=NF; i++) {
          if ($i=="A1")   cA1=i
          if ($i=="A2")   cA2=i
          if ($i=="AF1")  cAF1=i
          if ($i=="BETA") cBETA=i
        }
        print
        next
      }
      $1==CHR {
        # A1とA2 を交換し、アレル頻度を1から引いて、BETAの正負を変える
        tmp = $cA1
        $cA1 = $cA2
        $cA2 = tmp

        $cAF1 = 1 - $cAF1

        $cBETA = -$cBETA

        print
      }' \
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
      | sort -k1,1 \
      > $TEMP2

    zcat inputs/geneatlas/snps.imputed.chr${CHR}.csv.gz \
      | awk '
          {print $1, $2, $3, $4}
        ' \
      | sort -k1,1 \
      > $TEMP3

    {
      echo -e "SNP ALLELE iscores NBETA NSE PV Position A1 A2"
      join -1 1 -2 1 "$TEMP2" "$TEMP3"
    } \
      | tr ' ' '\t' \
      > $TEMP1
    ;;

  pan-ukb_*)
    # ここはもっとシンプルに書けるが、そうすると何故かエラーを吐かずに止まって解決できなかったので、動くがややこい書き方になっている、、
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
    #zcat $GWA_SUMSTA > $TEMP1
    zcat "$GWA_SUMSTA" \
      | awk -v CHR="$CHR" 'NR==1 || $1 ~ ("^" CHR ":")' \
      > "$TEMP1"
    ;;

  *)
    # tempファイル名から圧縮の有無を推測できないので解凍しておく
    ( [[ $GWA_SUMSTA == *.gz ]] && zcat $GWA_SUMSTA || cat $GWA_SUMSTA ) > $TEMP1
    ;;
esac

# ---- 最終的な出力を GWA_SUMSTA として置き換える ----
GWA_SUMSTA=$TEMP1


# 1: gwas sumstaをキレイにする
cp $GWA_SUMSTA aaaaaaaa.txt
python scripts/01_munge_sumsta_by_gwaslab.py \
  --sumsta $GWA_SUMSTA \
  --fmt $GWA_FMT \
  --rsid $LEAD_RS \
  --threads $THREAD \
  --out $GWA_SUMSTA_SUBSET \
  #--window $WINDOW \
  #> logs/1.$ANALYSIS.log 2>&1


## ---- N, N_CASE, TRAIT_TYPE ----
#case "$GWA_FMT" in
#  gwascatalog_hm|ssf)
#    eval "$(
#      python modules/get_sumsta_meta_gwascatalog.py \
#        --gwa-sumsta "$GWA_SUMSTA_BASE" \
#        --traits-table inputs/gwas_catalog-v1.0.3.1-studies-r2025-12-22.tsv
#    )"
#    ;;
#
#  fastgwa)
#    eval "$(
#      python modules/get_sumsta_meta_fastgwa.py \
#        --gwa-sumsta "$GWA_SUMSTA_BASE"
#    )"
#    ;;
#
#  geneatlas)
#    eval "$(
#      python modules/get_sumsta_meta_geneatlas.py \
#        --gwa-sumsta "$GWA_SUMSTA_BASE" \
#        --traits-table inputs/Traits_Table_GeneATLAS.csv \
#        --default-n 452264
#    )"
#    ;;
#
#  vcf)
#    # 実はopengwas
#    eval "$(
#      python modules/get_sumsta_meta_opengwas.py \
#        --gwa-sumsta "$GWA_SUMSTA_BASE" \
#        --traits-table inputs/opengwas.io_2026-01-05T07-03-28.csv
#    )"
#    ;;
#
#  gwasatlas)
#    eval "$(
#      python modules/get_sumsta_meta_manual.py \
#        --gwa-sumsta "$GWA_SUMSTA_BASE" \
#        --traits-table inputs/trait_info_manual_gwasatlas.txt
#    )"
#    ;;
#
#  pheweb)
#    # 実はfinngen
#    eval "$(
#      python modules/get_sumsta_meta_finngen.py \
#        --gwa-sumsta "$GWA_SUMSTA_BASE" \
#        --traits-table inputs/endpoints_finngen_r9.tsv
#    )"
#    ;;
#
#  pan-ukb_*)
#    eval "$(
#      python modules/get_sumsta_meta_panukbb.py \
#        --gwa-sumsta "$GWA_SUMSTA_BASE" \
#        --traits-table inputs/Pan-UK_Biobank_phenotype_manifest.txt
#    )"
#    ;;
#
#  NealeLab)
#    eval "$(
#      python modules/get_sumsta_meta_nealelab.py \
#        --gwa-sumsta "$GWA_SUMSTA_BASE" \
#        --traits-table inputs/phenotypes.both_sexes.v2_nealelab.tsv.bgz
#    )"
#    ;;
#
#  *)
#    eval "$(
#      python modules/get_sumsta_meta_manual.py \
#        --gwa-sumsta "$GWA_SUMSTA_BASE" \
#        --traits-table inputs/trait_info_manual_others.txt
#    )"
#    ;;
#esac
#
#echo $GWA_SUMSTA_BASE $TRAIT_TYPE $N $N_CASE $N_CONTROL
#
#if [[ "$TRAIT_TYPE" == *check* ]] || \
#   [[ "$N" == *check* ]] || \
#   [[ "$N_CASE" == *check* ]] || \
#   [[ "$N_CONTROL" == *check* ]]; then
#  eval "$(
#    python modules/get_sumsta_meta_manual.py \
#      --gwa-sumsta "$GWA_SUMSTA_BASE" \
#      --traits-table inputs/trait_info_manual.txt
#  )"
#fi
#
#echo $GWA_SUMSTA_BASE $TRAIT_TYPE $N $N_CASE $N_CONTROL



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
#
#
## 4: gwaslabでgtfとプロットする
## to do: 中身が01とまろかぶりなので、その部分をsumstaの共通関数とする
#GTF=inputs/gm24_gencodev47.chr.scaffold.444141isoforms.geneid_corrected.geneid_refined.add_txType_geneType_geneRow.sort.noChr.gtf
#
#python scripts/04_regional_plot_by_gwaslab.py \
#  --sumsta $GWA_SUMSTA \
#  --fmt $GWA_FMT \
#  --rsid $LEAD_RS \
#  --gtf $GTF \
#  --pop $LD_POP \
#  --window $WINDOW \
#  --threads $THREAD \
#  --out results/04/$LEAD_RS.$ENSG.$(basename $GWA_SUMSTA_ORIG).$LD_POP.gwaslab.png
