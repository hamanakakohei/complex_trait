#!/usr/bin/env bash
# 
# 注意点：
# 全数もQTL数もcausaldb側とjoinするので、VAR_TYPEはcausaldb側のみ絞っている
#
# 入力：
# rnaseq430.cis_independent_qtl.txt.gz: tensorqtl indepの結果
# gene_category_map.txt: gene_id, structural_category列の2列からなる
# credible_set.txt.gz: CAUSALdb提供のファイル
# meta.txt: CAUSALdb提供のファイル
# bin_scheme.txt
# rsid_variant__qtl.txt.gz: rsid、variant_id列の2列からなる、QTLで解析対象（e.g. TSS+-1Mb）のバリアントリスト、事前に作っている
#
# 出力：
# 4:
# qtl_in_filtered_causaldb.xxx.meta_info.txt: CAUSALdbの高PIP群（指定可）に落ちたQTLの情報（CAUSALdbのmeta.txt貼り付け）
# url_trait.txt: url, trait列の2列からなる、上のファイル中のGWAS sumstaをダウンロードするために必要な情報
#
# 5:
# フォレストプロット

set -euo pipefail
eval "$(conda shell.bash hook)"

conda activate misc_20250301


QTL_TYPE="e"
QTL_INDEP_SUMSTA=inputs/rnaseq430.cis_independent_qtl.txt.gz
QTL_INDEP=""
GENE_TYPE_MAP=inputs/gene_category_map.txt
CAUSALDB_CS=inputs/credible_set.txt.gz
CAUSALDB_INDEP=""
CAUSALDB_META=inputs/meta.txt
METHOD=""
VAR_TYPE=""
BIN_SCHEME=inputs/bin_scheme.txt
QTL_ALL_VAR=inputs/rsid_variant__qtl.txt.gz
RSID_VARID_MAP=inputs/rsid_variant__qtl.txt.gz

while [[ $# -gt 0 ]]; do
  case "$1" in
    --qtl_type)       QTL_TYPE="$2";       shift 2;;
    --qtl_indep)      QTL_INDEP="$2";      shift 2;;
    --method)         METHOD="$2";         shift 2;;
    --var_type)       VAR_TYPE="$2";       shift 2;;
    --causaldb_indep) CAUSALDB_INDEP="$2"; shift 2;;
    *) echo "Unknown option: $1";  exit 1;;
  esac
done


ANALYSIS_TYPE=$QTL_TYPE.$QTL_INDEP.$METHOD.$VAR_TYPE.$CAUSALDB_INDEP
GENE_TYPES=(known new antisense genic_intron intergenic)


# 1: causaldbの表を指定した条件でフィルターしつつ、バリアントーmaxPIPの表を作る
OPTS=""
if [ "$VAR_TYPE" = "snp" ];       then OPTS="$OPTS --snp_only"; echo "snp only"; fi
if [ "$CAUSALDB_INDEP" = "PrimGwas" ]; then OPTS="$OPTS --primary_only"; echo "CAUSALdb primary only"; fi

scripts/01.py \
  --causaldb $CAUSALDB_CS \
  --method $METHOD \
  --cs95_only \
  $OPTS \
  --prefix results/$ANALYSIS_TYPE/01/causaldb.filtered
  
  
# 2-1 & 2-2
for GENE_TYPE in ${GENE_TYPES[@]}; do
  # 2-1: qtl側を指定した条件（gene classとprimary signalかどうか）でフィルターする
  OPTS=""
  if [ "$QTL_INDEP" = "PrimQTL" ]; then OPTS="$OPTS --primary_only"; echo "QTL primary only"; fi
  
  scripts/02-1.py \
    --qtl_indep $QTL_INDEP_SUMSTA \
    --gene_class_map $GENE_TYPE_MAP \
    --gene_class $GENE_TYPE \
    --rsid_varid_map $RSID_VARID_MAP \
    $OPTS \
    --out results/$ANALYSIS_TYPE/02/qtl.filtered.$GENE_TYPE.txt


  # 2-2: PIPビンごとのQTL or notなバリアント数の分割表作る
  # 全数は１でフィルターしたcausaldb
  # このQTLファイルを事前にフィルターしてしまったせいでわかりにくくなったかも、、、
  # ついでに、QTLで、causaldb内にあるもの（pipでフィルターしない、0.0とかもある）リストを出力している
  scripts/02-2.py \
    --qtl results/$ANALYSIS_TYPE/02/qtl.filtered.$GENE_TYPE.txt \
    --causaldb_pip results/$ANALYSIS_TYPE/01/causaldb.filtered.rsid156_vs_maxPIP.txt.gz \
    --bin_table $BIN_SCHEME \
    --qtl_analysis_variants $QTL_ALL_VAR \
    --prefix results/$ANALYSIS_TYPE/02/qtl_in_filtered_causaldb.$GENE_TYPE
done


# 2-3：遺伝子クラスを全てまとめて、"AnyGene"とする
# 5でプロット用に使う
OUT=results/$ANALYSIS_TYPE/02/qtl_in_filtered_causaldb.AnyGene.bin_count.txt
echo -e "bin\tn_in_bin\tn_out_bin\tdata" > $OUT

for GENE_TYPE in ${GENE_TYPES[@]}; do
  IN=results/$ANALYSIS_TYPE/02/qtl_in_filtered_causaldb.$GENE_TYPE.bin_count.txt
  awk -v GENE_TYPE=${GENE_TYPE} 'NR>1{print $0"\t"GENE_TYPE}' $IN >> $OUT
done


# 3：03.pyにMeSH IDでcausaldbをフィルターするやり方を作っている途中


# 4：new遺伝子のQTLにcausaldbのメタ情報を付けたファイルを作る
# ここ２で作った、QTLとcausaldbのかぶりリストを利用しているが、わかりにくいかもしれない
# QTLとcausaldbファイルを与えて被りリストを内部で作り直したほうがわかりやすいか
scripts/04.py \
  --causaldb_cs results/$ANALYSIS_TYPE/01/causaldb.filtered.whole.txt.gz \
  --causaldb_fm_method pip \
  --causaldb_pip_thr 0.1 \
  --target_rsid_file results/$ANALYSIS_TYPE/02/qtl_in_filtered_causaldb.new.txt.gz \
  --causaldb_meta $CAUSALDB_META \
  --qtl_meta $QTL_INDEP_SUMSTA \
  --rsid_varid_map $RSID_VARID_MAP \
  --out1 results/$ANALYSIS_TYPE/04/qtl_in_filtered_causaldb.new.meta_info.txt \
  --out2 results/$ANALYSIS_TYPE/04/url_trait.txt
  ## 01のフィルター済みCSファイルを使わない場合は以下：
  #--causaldb_fm_method $METHOD \
  #--causaldb_cs $CAUSALDB_CS \



# 5：フォレストプロット（2-2で作ったAnyGene.bin_countファイルを基にする）
conda activate misc_r

scripts/05.R \
  --contingency_table results/$ANALYSIS_TYPE/02/qtl_in_filtered_causaldb.AnyGene.bin_count.txt \
  --effect_col n_in_bin \
  --noneffect_col n_out_bin \
  --ref_bin b_0.0_0.01 \
  --bin_order b_0.0_0.01 b_0.01_0.1 b_0.1_0.5 b_0.5_1.0 \
  --data_order known new antisense genic_intron intergenic \
  --out results/$ANALYSIS_TYPE/05/AnyGene.png
