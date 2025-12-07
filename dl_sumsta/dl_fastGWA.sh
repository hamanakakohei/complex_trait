#!/usr/bin/env bash
#
# fastGWA sumstaをダウンロードする.
# URLが更新されて候補となるのが複数あるので、与えたURLからファイル名を取り出して全部試す
# なのでset -e は外している（失敗時に次へ進めないため）
#
# Usage:
#   dl_fastGWA.sh <URL> <OUT_DIR> <THREADS>
#
# Arguments:
#   URL        Primary URL to try first (others are auto-generated).
#   OUT_DIR    Directory where the downloaded file will be saved.
#   THREADS    Number of threads passed to aria2c (-x and -s).
#
# Notes:
#   - Output directory is created automatically if missing.
#   - このスクリプトの出力はタブ区切りで:
#       "0 <original_URL> <filename>" （成功時）
#       "1 <original_URL>" （失敗時）
#----------------------------------------------
set -u
set -o pipefail

URL="$1"
OUT_DIR="$2"
THREAD="$3"

mkdir -p "$OUT_DIR"


# URLの候補
f=$(basename "$URL")
URL2="https://yanglab.westlake.edu.cn/data/fastgwa_data/UKB/$f"
URL3="https://yanglab.westlake.edu.cn/data/fastgwa_data/UKBbin/$f"
URL4="https://yanglab.westlake.edu.cn/data/fastgwa_data/WES/$f"
URLS=("$URL" "$URL2" "$URL3" "$URL4")


SUCCESS=1

for U in "${URLS[@]}"; do
    if aria2c \
        --auto-file-renaming=false \
        --continue=true \
        -q \
        -x "$THREAD" \
        -s "$THREAD" \
        -j 1 \
        -d "$OUT_DIR" \
        "$U"; then

        SUCCESS=0
        break
    fi
done


# 全URLで失敗なら1を返す
if [[ $SUCCESS -eq 0 ]]; then
    echo -e "0\t$URL\t$f"
    exit 0
else
    echo -e "1\t$URL"
    exit 1
fi
