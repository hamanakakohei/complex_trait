#!/bin/bash
#SBATCH --job-name=qtl2gwas_coloc
#SBATCH --partition=epyc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate coloc3

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1


INPUT=$1
EQTL_ASSOC_PREFIX=$2
EQTL_GENO_PREFIX=$3
EQTL_PHEN=$4
EQTL_COV=$5

IFS=$'\t' read -r RSID VARID CHR EGENE F FMT GWAS_LD_POP URL TRAIT MESH_TERM MESH_ID GWAS_N GWAS_N_CASE GWAS_N_CONTROL GWAS_TYPE \
        MY_GWAS_N MY_GWAS_N_CASE MY_GWAS_N_CONTROL \
        < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $INPUT)

GWAS_SUMSTA=results/01/${F}.chr${CHR}.txt.gz
GWAS_N=${GWAS_N%.0}
GWAS_N_CASE=${GWAS_N_CASE%.0}
GWAS_N_CONTROL=${GWAS_N_CONTROL%.0}
MY_GWAS_N_CASE=${MY_GWAS_N_CASE%.0}

WINDOW=2000000

mkdir -p logs


GWAS_N_CASE_ARG=""
if [ "$GWAS_TYPE" = "cc" ]; then
  if [ "$GWAS_N_CASE" = "NA" ]; then
    GWAS_N_CASE_ARG="--gwas_n_case $MY_GWAS_N_CASE"
  else
    GWAS_N_CASE_ARG="--gwas_n_case $GWAS_N_CASE"
  fi
fi

if [[ "$GWAS_LD_POP" == "UKB" ]]; then
  GWAS_LD_POP="EUR"
fi

echo $VARID
echo $EGENE
echo $EQTL_ASSOC_PREFIX
echo $EQTL_GENO_PREFIX
echo $EQTL_PHEN
echo $EQTL_COV
echo "window: $WINDOW"
echo $GWAS_SUMSTA
echo $GWAS_TYPE
echo $GWAS_LD_POP
echo "gwas_n: $GWAS_N"
echo "gwas_n_case: $GWAS_N_CASE"
echo "gwas_n_control: $GWAS_N_CONTROL"


Rscript scripts/02_qtl_gwas_coloc.R \
        --variant_id $VARID \
        --egene $EGENE \
        --eqtl_assoc_prefix $EQTL_ASSOC_PREFIX \
        --eqtl_geno_prefix $EQTL_GENO_PREFIX \
        --eqtl_phen $EQTL_PHEN \
        --eqtl_cov $EQTL_COV \
        --window $WINDOW \
        --gwas_sumsta $GWAS_SUMSTA \
        --gwas_type $GWAS_TYPE \
        --gwas_ld_pop $GWAS_LD_POP \
        --gwas_n $GWAS_N \
        $GWAS_N_CASE_ARG \
        --out_prefix results/part2/02/${SLURM_ARRAY_TASK_ID}.${WINDOW}
        #--only_plot \


#SUMMARY_FILE=$1
#QTL_SUMSTA_PREFIX=$2
#
#THREADS=${SLURM_CPUS_PER_TASK:-1}
#
#
## 準備
#IFS=$'\t' read TRAIT_TYPE LEAD_RS CHR ENSG GWA_SUMSTA GWA_FMT LDPOP \
#  <<< "$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{OFS="\t"}NR==n {print $1, $2, $3, $4, $5, $6, $7}' "$SUMMARY_FILE")"
#
#if [[ "$LDPOP" == "UKB" ]]; then
#    LDPOP="EUR"
#fi
#
#QTL_SUMSTA=${QTL_SUMSTA_PREFIX}${CHR}.parquet
#
#
#echo "LEAD_RS: $LEAD_RS"
#echo "CHR: $CHR"
#echo "ENSG: $ENSG"
#echo "GWA_SUMSTA: $GWA_SUMSTA"
#echo "GWA_FMT: $GWA_FMT"
#echo "LDPOP: $LDPOP"
#echo "QTL_SUMSTA: $QTL_SUMSTA"
#echo "THREADS: $THREADS"
#
#
## plot_windowは本当はgwas sumstaでリードを中心に限定するのみだが、
## locuscompareでqtl側と内部結合するので結局プロットウィンドウとなる
## ここで指定したLEADRSがデータ内になければstop()する？
## ENSGはQTL SUMSTAのphenotypeを限定するため
#bash pipeline/pipeline_qtl_gwas_coloc.sh \
#  --lead_rs $LEAD_RS \
#  --chr $CHR \
#  --gwa_sumsta inputs/gwas_sumsta/$GWA_SUMSTA \
#  --gwa_fmt $GWA_FMT \
#  --qtl_sumsta $QTL_SUMSTA \
#  --ensg $ENSG \
#  --plot_window 300 \
#  --ld_pop $LDPOP \
#  --thread $THREADS \
#  --out results/03/$LEAD_RS.$ENSG.$(basename $GWA_SUMSTA).$LDPOP.png
