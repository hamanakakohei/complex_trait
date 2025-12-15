#!/bin/bash
#SBATCH --job-name=01_locuscompare
#SBATCH --partition=epyc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

set -euo pipefail

SUMMARY_FILE=$1
QTL_SUMSTA_PREFIX=$2

THREADS=${SLURM_CPUS_PER_TASK:-1}


# 準備
IFS=$'\t' read RS CHR ENSG GWA_SUMSTA GWA_FMT LDPOP \
  <<< "$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{OFS="\t"}NR==(n+1) {print $1, $2, $3, $4, $5, $6}' "$SUMMARY_FILE")"

if [[ "$LDPOP" == "UKB" ]]; then
    LDPOP="EUR"
fi

QTL_SUMSTA=${QTL_SUMSTA_PREFIX}${CHR}.parquet


echo "RS: $RS"
echo "CHR: $CHR"
echo "ENSG: $ENSG"
echo "GWA_SUMSTA: $GWA_SUMSTA"
echo "GWA_FMT: $GWA_FMT"
echo "LDPOP: $LDPOP"
echo "QTL_SUMSTA: $QTL_SUMSTA"
echo "THREADS: $THREADS"


# plot_windowは本当はgwas sumstaでリードを中心に限定するのみだが、
# locuscompareでqtl側と内部結合するので結局プロットウィンドウとなる
# ここで指定したLEADRSがデータ内になければstop()する？
# ENSGはQTL SUMSTAのphenotypeを限定するため
bash pipeline/pipeline.sh \
  --lead_rs $RS \
  --chr $CHR \
  --ensg $ENSG \
  --gwa_sumsta inputs/gwas_sumsta/$GWA_SUMSTA \
  --gwa_fmt $GWA_FMT \
  --qtl_sumsta $QTL_SUMSTA \
  --plot_window 300 \
  --ld_pop $LDPOP \
  --thread $THREADS \
  --out results/03/$RS.$ENSG.$(basename $GWA_SUMSTA).$LDPOP.png
