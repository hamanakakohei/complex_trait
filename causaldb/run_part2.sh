#!/usr/bin/env bash
#
# 01
#   入力
#     gcst_dir_patch.txt:
#       タブ区切りで置換前、置換後を並べたもの
#       02で失敗してこのパッチを更新していった
#       gwas catalogのURLが更新前のままなのを直したり、pubmedのをgwas catalogのurlに変えたり
#
# WARNING
#   02でDL出来なかったsumsta
#     http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630
#     http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008916
#     bloodcellgenetics
#
# 03
#   入力
#     url_trait_DLedFile.txt
#       - url, trait, f: この3列は02のログの標準出力のラスト1行をcatして作る
#       - fmt: 各sumstaをgwaslabで読むときのフォーマットを手入力する、適切なものがなければformatbook.jsonを加筆する
#
#   出力
#     locuscompare_input.txt
#       - rsid: eQTL
#       - chr: eQTLの染色体番号で、gwas sumstaをこれでフィルターして読む
#       - phenotype_id: eGene
#       - f: gwas sumstaのファイル名
#       - fmt: それをgwaslabで読む時に指定するフォーマット名（geneatlas, fastgwa, etc.
#       - popu: locuscompareRでregional plotするときに使うLD用POP
#       - url: CAUSALdb提供でsumstaをDLするときに使ったもの
#       - trait: CAUSALdb提供でsumstaをDLするときに使ったもの
set -euo pipefail
eval "$(conda shell.bash hook)"

conda activate misc_20250301


while [[ $# -gt 0 ]]; do
  case "$1" in
    --analysis_type) ANALYSIS_TYPE="$2";      shift 2;;
    *) echo "Unknown option: $1";  exit 1;;
  esac
done


# 1: URLなどを手直しした方が早い場合はここで直す
sed -f <(
    awk -F'\t' '{printf "s|%s|%s|g\n", $1, $2}' inputs/gcst_dir_patch.txt
  ) results/$ANALYSIS_TYPE/part1/04/url_trait.txt \
  > results/$ANALYSIS_TYPE/part2/01/url_trait.patched.txt

sed -f <(
    awk -F'\t' '{printf "s|%s|%s|g\n", $1, $2}' inputs/gcst_dir_patch.txt
  ) results/$ANALYSIS_TYPE/part1/04/qtl_in_filtered_causaldb.new.meta_info.txt \
  > results/$ANALYSIS_TYPE/part2/01/qtl_in_filtered_causaldb.new.meta_info.patched.txt


# 2: gwas sumstaをダウンロードする
N_SUMSTA=$(wc -l < results/$ANALYSIS_TYPE/part2/01/url_trait.patched.txt)

#sbatch --array=1,10,53,66,68,72,133,138,147,148,163%6 \
sbatch --array=1-"$N_SUMSTA"%5 \
  slurm/part2/02_dl_sumsta_array.slurm \
  results/$ANALYSIS_TYPE/part2/01/url_trait.patched.txt \
  results/$ANALYSIS_TYPE/part2/02/


# 3: locuscompareの入力とするファイルを作る
scripts/part2/03.py \
  --in1 results/$ANALYSIS_TYPE/part2/01/qtl_in_filtered_causaldb.new.meta_info.patched.txt \
  --in2 results/$ANALYSIS_TYPE/part2/02/url_trait_DLedFile.txt \
  --out results/$ANALYSIS_TYPE/part2/03/locuscompare_input.txt
