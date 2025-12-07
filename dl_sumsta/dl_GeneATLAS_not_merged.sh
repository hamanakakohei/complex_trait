#!/usr/bin/env bash
#
# GeneATLAS sumstaをダウンロードし。全染色体をマージせずにそのまま
# マージするとファイルが巨大かつ各行の染色体が不明になり後でフィルターできないので
# 
#
# Usage:
#   dl_GeneATLAS_not_merged.sh <TRAIT> <MANIFEST> <OUT_DIR>
#
# Arguments:
#   TRAIT        例：https://gwas.mrcieu.ac.uk/files/ukb-b-6027/ukb-b-6027.vcf.gz
#   MANIFEST   http://geneatlas.roslin.ed.ac.uk/traits-table/からdlできるTraits_Table_GeneATLAS.csvのこと
#   OUT_DIR    Directory where the downloaded file will be saved.
#
# Notes:
#   - Output directory is created automatically if missing.
#   - このスクリプトの出力はタブ区切りで:
#       "0 <trait> <file_prefix>" （成功時）
#       "1 <trait>" （失敗時）
#   - この<file_prefix>は、${OUT_DIR}/imputed.allWhites.${KEY}.chrみたいな感じで、後で染色体番号を指定して.csv.gzを付けて、ファイルパスが完成する
#----------------------------------------------
set -euo pipefail

TRAIT="$1"
MANIFEST="$2"
OUT_DIR="$3"

mkdir -p "$OUT_DIR"


# ダウンロード用にTRAIT --> KEY --> URLとしないといけない
# MANIFESTの2列目のDescription列の値がTRAITと一致する行の、1列目のkey列の値を取り出す
# 最後にtail -n1しているのは、WeightとBMIが2行（2 key）あるため、、（mlr --csv cut -f Description $MANIFEST | sort | uniq -cで確認できる）
KEY=$(mlr --csv filter "\$Description==\"$TRAIT\"" then cut -f key $MANIFEST | tail -n1)

url_for_chr() {
    local chr="$1"
    if [[ "$chr" == "X" ]]; then
        echo "http://static.geneatlas.roslin.ed.ac.uk/gwas/allWhites/imputed/data.chromX/base/imputed.allWhites.combined.${KEY}.chrX.csv.gz"
    else
        echo "http://static.geneatlas.roslin.ed.ac.uk/gwas/allWhites/imputed/data.copy/imputed.allWhites.${KEY}.chr${chr}.csv.gz"
        #echo "http://static.geneatlas.roslin.ed.ac.uk/gwas/allWhites/genotyped/data/genotyped.allWhites.${KEY}.chr${chr}.csv.gz"
    fi
}


# 各染色体のファイルをダウンロード
for CHR in {1..22} X; do
    URL=$(url_for_chr "$CHR")
    FILE="${OUT_DIR}/imputed.allWhites.${KEY}.chr${CHR}.csv.gz"

    # ダウンロード
    wget -c -q -O "$FILE" "$URL" || {
        echo -e "1\t$TRAIT"
        exit 1
    }
done


# prefix を抽出: chr の直前まで
prefix="${OUT_DIR}/imputed.allWhites.${KEY}.chr"
echo -e "0\t$TRAIT\t$prefix"
