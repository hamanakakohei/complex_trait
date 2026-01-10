#!/bin/bash
#SBATCH --job-name=qtl2qtl_coloc
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#
# 注:
# regionalプロットを描きたいだけの時は--only_plotを付ける
# genotypeをcovariateでregress outしてからin-sample LDを計算したい時はCOV_ADJ_LDを指定する
# coloc対象領域をWINDOWで指定できる
#
# to do:
# これらをすべてこのスクリプトに対する引数で指定出来るようにする
set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate coloc2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

INPUT_SUMMARY=$1
ASSOC_PREFIX=$2
GENO_PREFIX=$3
PHEN=$4
COV=$5

WINDOW=500000
COV_ADJ_LD=true

read -r PHEN1 LEAD RANK PHEN2 OTHERS < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $INPUT_SUMMARY)


COV_ADJ_LD_FLAG=""
if [ "$COV_ADJ_LD" = true ]; then
  COV_ADJ_LD_FLAG="--adjust_ld_by_cov"
fi

Rscript scripts/01_qtl_qtl_coloc.R \
        --variant_id $LEAD \
        --phen_id1 $PHEN1 \
        --phen_id2 $PHEN2 \
        --assoc_prefix $ASSOC_PREFIX \
        --geno_prefix $GENO_PREFIX \
        --phen $PHEN \
        --cov $COV \
        --window $WINDOW \
        $COV_ADJ_LD_FLAG \
        --out_prefix results/01/${SLURM_ARRAY_TASK_ID}.window_${WINDOW}.covAdj_$COV_ADJ_LD
        #--only_plot 
