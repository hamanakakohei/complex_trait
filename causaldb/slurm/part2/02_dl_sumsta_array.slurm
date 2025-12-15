#!/bin/bash
#SBATCH --job-name=dl_sumsta
#SBATCH --partition=epyc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=6G
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

set -euo pipefail


URL_TRAIT=$1
OUT_DIR=$2

IFS=$'\t' read -r URL TRAIT < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$URL_TRAIT")


bash pipeline/part2/02_dl_sumsta.sh \
  --url $URL \
  --trait "$TRAIT" \
  --out_dir $OUT_DIR
