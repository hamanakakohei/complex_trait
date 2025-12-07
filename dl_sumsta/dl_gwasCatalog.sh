#!/usr/bin/env bash
#
# Downloads harmonised GWAS summary statistics (.h.tsv.gz) from a given
# (ftp or https) GWAS catalog URL. Automatically converts ftp:// to
# https:// and searches the `/harmonised/` directory for GRCh38 files.
#
# Usage:
#   dl_gwas_harmonised.sh <URL> <OUTDIR> [THREADS]
#
# Arguments:
#   URL         Base URL (ftp:// or https://)
#   OUTDIR      Output directory for downloaded files
#   THREADS     Number of threads for aria2c (default: 5)
#
# Output:
#   0  URL  FILENAME   ... success
#   1  URL             ... failure
#
# Notes:
#   * The script retrieves the file list from <URL>/harmonised/
#   * Searches for files ending with ".h.tsv.gz" (GRCh38 harmonised)

set -euo pipefail

URL=$1
OUTDIR=$2
THREAD=${3:-5}

mkdir -p "$OUTDIR"


# --- ftp:// → https:// に置換 ---
URL2="${URL/ftp:/https:}"
HARM_URL="${URL2}/harmonised"


# --- harmonised ディレクトリを-listing ---
page=$(wget -q -O - "${HARM_URL}/" || true)

if [[ -z "$page" ]]; then
    echo -e "1\t$URL"
    exit 1
fi


# --- ファイル一覧抽出 ---
file_list=$(echo "$page" | grep -oE 'href="[^"]+"' | cut -d'"' -f2 || true)

if [[ -z "$file_list" ]]; then
    echo -e "1\t$URL"
    exit 1
fi


# --- GRCh38 .h.tsv.gz を抽出 ---
target=$(echo "$file_list" | grep -E '\.h\.tsv\.gz$' || true)

if [[ -z "$target" ]]; then
    echo -e "1\t$URL"
    exit 1
fi


# --- ダウンロード ---
for f in $target; do
    full="${HARM_URL}/${f}"

    # aria2c による並列ダウンロード
    if aria2c --auto-file-renaming=false --continue=true -q -s "$THREAD" -x "$THREAD" -j 1 -d "$OUTDIR" "$full"; then
        echo -e "0\t$URL\t$f"
        exit 0
    else
        echo -e "1\t$URL"
        exit 1
    fi
done
