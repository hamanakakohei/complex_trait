#!/usr/bin/env bash
#
# GeneATLAS sumstaをダウンロードし。全染色体分をマージする
# 詳しくは、dl_GeneATLAS_not_merged.shと同じ
#
# Usage:
#   dl_GeneATLAS_merged.sh <TRAIT> <MANIFEST> <OUT_DIR>
#
# Arguments:
#   TRAIT      例："Ease of skin tanning"
#   MANIFEST   http://geneatlas.roslin.ed.ac.uk/traits-table/からdlできるTraits_Table_GeneATLAS.csvのこと
#   OUT_DIR    Directory where the downloaded file will be saved.
#
# Notes:
#   - Output directory is created automatically if missing.
#   - このスクリプトの出力はタブ区切りで:
#       "0 <trait> <file_name>" （全ファイル成功時）
#       "1 <trait>" （一つでも失敗時）
#   - この<file_name>は、${KEY}.merged.csv.gz
#----------------------------------------------
set -euo pipefail

TRAIT="$1"
MANIFEST="$2"
OUT_DIR="$3"

mkdir -p "$OUT_DIR"


# ダウンロード用にTRAIT --> KEY --> URLとしないといけない
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


# 各染色体のファイルをダウンロードしつつ、結合していく
MERGED="${OUT_DIR}/${KEY}.merged.csv.gz"
TMPDIR=$(mktemp -d)

first=1

for CHR in {1..22} X; do
    URL=$(url_for_chr "$CHR")
    FILE="${TMPDIR}/chr${CHR}.csv.gz"

    # ダウンロード
    wget -c -q -O "$FILE" "$URL" || {
        echo -e "1\t$TRAIT"
        exit 1
    }

    # 結合処理（ヘッダーは1回だけ）
    if [[ $first -eq 1 ]]; then
        zcat "$FILE" > "${TMPDIR}/merged.csv"
        first=0
    else
        zcat "$FILE" | tail -n+2 >> "${TMPDIR}/merged.csv"
    fi
done


# 保存する
gzip -c "${TMPDIR}/merged.csv" > "$MERGED"

f=$(basename $MERGED)
echo -e "0\t$TRAIT\t$f"
rm -r "$TMPDIR"
