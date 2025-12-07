#!/usr/bin/env bash
#
# openGWAS sumstaをダウンロードする.
# トークンが要るので、.bash_profileなどにexport OPENGWAS_TOKEN=xxxを書いて、事前にsource .bash_profileとかしておく
# 今は、CAUSALdbの示す間違ったURLからtrait? IDを抜き出して、APIで真のURLを得ているのだが、本当は引数でIDを直に与えた方がシンプル
#
# Usage:
#   dl_openGWAS.sh <URL> <OUT_DIR> <THREADS>
#
# Arguments:
#   URL        例：https://gwas.mrcieu.ac.uk/files/ukb-b-6027/ukb-b-6027.vcf.gz
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
