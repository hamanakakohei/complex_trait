#!/usr/bin/env bash
#
# openGWAS sumstaをダウンロードする.
# URLが更新されて候補となるのが複数あるので、与えたURLからファイル名を取り出して全部試す
#
# Usage:
#   download_fallback.sh <URL> <OUT_DIR> <THREADS>
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
set -euo pipefail

URL="$1"
OUT_DIR="$2"
THREAD="$3"

mkdir -p "$OUT_DIR"


# WARNING: URL の最後から2番目がID になる想定
# でAPIで真の一時的なURLを得る
ID=$(awk -F"/" '{print $(NF-1)}' <<< "$URL")
API="https://api.opengwas.io/api/gwasinfo/files?id=${ID}&id=string"

json=$(curl -sS -X 'POST' \
    $API \
    -H 'accept: application/json' \
    -H "Authorization: Bearer $OPENGWAS_TOKEN" \
    -d '')


# そこからVCFのURLを選んでダウンロード
vcf_url=$(jq -r ".[\"$ID\"][] | select(endswith(\".vcf.gz\"))" <<< "$json" | head -n 1)

if aria2c --auto-file-renaming=false --continue=true -q -x $THREAD -s $THREAD -j 1 -d $OUT_DIR $vcf_url; then
    f=$(basename $vcf_url)
    echo -e "0\t$URL\t$f"
    exit 0
else
    echo -e "1\t$URL"
    exit 1
fi
