#!/usr/bin/env bash
#
# Description: 
#   traitとNeale Labマニフェストファイルを与えて、
#   そのマニフェストからwgetコマンドを抜き出して実行して、OUT_DIRへmvする
#
# Usage:
#   ./dl_nealeUKBB_from_traits.sh <TRAIT> <MANIFEST> <OUT_DIR>
#
# Arguments:
#   TRAIT       Primary URL to try first (others are auto-generated).
#   OUT_DIR     Directory where the downloaded file will be saved.
#   MANIFEST    ukb31063_ldsc_sumstat_manifest_aws_sep2022.tsv
# 
# 出力:
#   0 <TRAIT> <filename>（成功時）
#   1 <TRAIT>（失敗時）
#
# Notes / Warnings:
#   This script relies on the manifest format:
#       col2: description
#       col4: is_primary_gwas (TRUE/FALSE)
#       col10: wget command
#       col12: output filename
#   Only primary GWAS entries (is_primary_gwas=TRUE) are used.
# ------------------------------
set -euo pipefail

TRAIT="$1"          
MANIFEST="$2"       
OUT_DIR="$3"

mkdir -p $OUT_DIR

awk -F"\t" -v trait="$TRAIT" -v out_dir="$OUT_DIR" '
BEGIN { OFS="\t" }
{
    if ($2 == trait && $4 == "TRUE") {

        cmd = $10" -c -q"
        filename = $12

        if (cmd == "NA" || cmd == "" || cmd !~ /^wget /) {
            # wget コマンドではない場合
            print "1", trait
            next
        }

        # 実行
        rc = system(cmd)

        if (rc != 0) {
            print "1", trait
            next
        }

        system("mv "filename" "out_dir"/")
        print "0", trait, filename
    }
}
' "$MANIFEST"
