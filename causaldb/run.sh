#!/usr/bin/env bash

eval "$(conda shell.bash hook)"
conda activate misc_20250301


# 1: causaldbの表をフィルターしつつ、バリアントーmaxPIPの表を作る
CAUSALDB=/path/to/causaldb/v2.1/credible_set.txt.gz
METHODS=(
  abf
  finemap
  paintor
  caviarbf
  susie
  polyfun_finemap
  polyfun_susie
)

for METHOD in ${METHODS[@]}; do
  scripts/01.py \
    --causaldb $CAUSALDB \
    --method $METHOD \
    --snp_only \
    --cs95_only \
    --primary_only \
    --prefix results/01/$METHOD.snp

  scripts/01.py \
    --causaldb $CAUSALDB \
    --method $METHOD \
    --cs95_only \
    --primary_only \
    --prefix results/01/$METHOD.all
done


# 2: PIPビンごとのQTL or notなバリアント数の分割表作る
VAR_TYPES=(snp all)
GENE_TYPES=(known new antisense genic_intron intergenic)
BIN_SCHEME=inputs/bin_scheme.txt
QTL_ALL_VAR=inputs/rsid_variant__qtl.txt.gz

for METHOD in ${METHODS[@]}; do
  for VAR_TYPE in ${VAR_TYPES[@]}; do
    for GENE_TYPE in ${GENE_TYPES[@]}; do
      QTL=inputs/e.${VAR_TYPE}.${GENE_TYPE}.txt

      scripts/02.py \
        --qtl $QTL \
        --causaldb_pip results/01/$METHOD.${VAR_TYPE}.rsid_156.maxPIP.txt.gz \
        --bin_table $BIN_SCHEME \
        --qtl_analysis_variants $QTL_ALL_VAR \
        --prefix results/02/$METHOD.e.${VAR_TYPE}.${GENE_TYPE}
    done
  done
done


# 2-2：遺伝子クラスを全てまとめて、"any"とする
for METHOD in ${METHODS[@]}; do
  for VAR_TYPE in ${VAR_TYPES[@]}; do
    OUT=results/02/$METHOD.e.${VAR_TYPE}.any.bin_count.txt
    echo -e "bin\tn_in_bin\tn_out_bin\tdata" > $OUT
    for GENE_TYPE in ${GENE_TYPES[@]}; do
      IN=results/02/$METHOD.e.${VAR_TYPE}.${GENE_TYPE}.bin_count.txt
      awk -v GENE_TYPE=${GENE_TYPE} 'NR>1{print $0"\t"GENE_TYPE}' $IN >> $OUT
    done
  done
done


# 3：フォレストプロット
conda activate misc_r

for METHOD in ${METHODS[@]}; do
  for VAR_TYPE in ${VAR_TYPES[@]}; do
    #for GENE_TYPE in ${GENE_TYPES[@]}; do
    for GENE_TYPE in "any"; do
      echo results/02/$METHOD.e.${VAR_TYPE}.${GENE_TYPE}.bin_count.txt
      scripts/03.R \
        --contingency_table results/02/$METHOD.e.${VAR_TYPE}.${GENE_TYPE}.bin_count.txt \
        --effect_col n_in_bin \
        --noneffect_col n_out_bin \
        --ref_bin b_0.0_0.01 \
        --bin_order b_0.0_0.01 b_0.01_0.1 b_0.1_0.5 b_0.5_1.0 \
        --data_order known new antisense genic_intron intergenic \
        --out results/03/$METHOD.e.${VAR_TYPE}.${GENE_TYPE}.png
    done
  done
done


# 4：04.pyにMeSH IDでcausaldbをフィルターするやり方を作っている途中
