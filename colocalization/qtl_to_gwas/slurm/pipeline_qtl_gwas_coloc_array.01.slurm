#!/bin/bash
#SBATCH --job-name=qtl2gwas_coloc
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
set -euo pipefail

SUMMARY_FILE=$1
QTL_SUMSTA_PREFIX=$2

THREADS=${SLURM_CPUS_PER_TASK:-1}


# 準備
IFS=$'\t' read LEAD_RS CHR GWA_SUMSTA GWA_FMT \
  <<< "$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR==n {print}' "$SUMMARY_FILE")"

#if [[ "$LDPOP" == "UKB" ]]; then
#    LDPOP="EUR"
#fi
#
#QTL_SUMSTA=${QTL_SUMSTA_PREFIX}${CHR}.parquet


echo "LEAD_RS: $LEAD_RS"
echo "CHR: $CHR"
#echo "ENSG: $ENSG"
echo "GWA_SUMSTA: $GWA_SUMSTA"
echo "GWA_FMT: $GWA_FMT"
#echo "LDPOP: $LDPOP"
#echo "QTL_SUMSTA: $QTL_SUMSTA"
echo "THREADS: $THREADS"


# plot_windowは本当はgwas sumstaでリードを中心に限定するのみだが、
# locuscompareでqtl側と内部結合するので結局プロットウィンドウとなる
# ここで指定したLEADRSがデータ内になければstop()する？
# ENSGはQTL SUMSTAのphenotypeを限定するため
bash pipeline/pipeline_qtl_gwas_coloc.sh \
  --lead_rs $LEAD_RS \
  --chr $CHR \
  --gwa_sumsta inputs/gwas_sumsta/$GWA_SUMSTA \
  --gwa_fmt $GWA_FMT \
  --thread $THREADS \
  #--qtl_sumsta $QTL_SUMSTA \
  #--ensg $ENSG \
  #--plot_window 300 \
  #--ld_pop $LDPOP \
  #--out results/03/$LEAD_RS.$ENSG.$(basename $GWA_SUMSTA).$LDPOP.png
