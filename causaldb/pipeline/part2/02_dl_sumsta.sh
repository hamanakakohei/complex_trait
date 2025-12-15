#!/usr/bin/env bash


# 「GCST」は他のdbも含んでることあるので一番最後にマッチを調べる
# bash_profileにOPENGWAS_TOKENを書いておく

# causaldbからだとtrait名でneale lab manifestファイルを検索するしかなさそう、、

# 以下はマニュアル
# GWASATLAS
# - trait\tsumstats
# - Comparative height size at age 10\tf.1697.0.0_res.EUR.sumstats.MACfilt.txt.gz
# - Male-specific factors - Hair/balding pattern: Pattern 1\t2395_1_logistic.EUR.sumstats.MACfilt.txt.gz
# - Pulse rate (automated reading)\tf.102.0.0_res.EUR.sumstats.MACfilt.txt.gz
# - Seen doctor (GP) for nerves, anxiety, tension or depression\tf.2090.0.0_logistic.EUR.sumstats.MACfilt.txt.gz
# - Sitting height\tf.20015.0.0_res.EUR.sumstats.MACfilt.txt.gz
# - Standing height\tf.50.0.0_res.EUR.sumstats.MACfilt.txt.gz
# - aria2c --auto-file-renaming=false --continue=true -q -s $THREAD -x $THREAD -j 1 -d $OUT_DIR https://atlas.ctglab.nl/ukb2_sumstats/f.1697.0.0_res.EUR.sumstats.MACfilt.txt.gz
# - aria2c --auto-file-renaming=false --continue=true -q -s $THREAD -x $THREAD -j 1 -d $OUT_DIR https://atlas.ctglab.nl/ukb2_sumstats/2395_1_logistic.EUR.sumstats.MACfilt.txt.gz
# - aria2c --auto-file-renaming=false --continue=true -q -s $THREAD -x $THREAD -j 1 -d $OUT_DIR https://atlas.ctglab.nl/ukb2_sumstats/f.102.0.0_res.EUR.sumstats.MACfilt.txt.gz
# - aria2c --auto-file-renaming=false --continue=true -q -s $THREAD -x $THREAD -j 1 -d $OUT_DIR https://atlas.ctglab.nl/ukb2_sumstats/f.2090.0.0_logistic.EUR.sumstats.MACfilt.txt.gz
# - aria2c --auto-file-renaming=false --continue=true -q -s $THREAD -x $THREAD -j 1 -d $OUT_DIR https://atlas.ctglab.nl/ukb2_sumstats/f.20015.0.0_res.EUR.sumstats.MACfilt.txt.gz
# - aria2c --auto-file-renaming=false --continue=true -q -s $THREAD -x $THREAD -j 1 -d $OUT_DIR https://atlas.ctglab.nl/ukb2_sumstats/f.50.0.0_res.EUR.sumstats.MACfilt.txt.gz
#
# JENGER (ブラウザーからのみ）
# - BBJ.AG.autosome.txt.gz
# - BBJ.NAP.autosome.txt.gz
# PUBMED
#
# bloodcellgenetics 

set -euo pipefail
eval "$(conda shell.bash hook)"

conda activate misc_20250301


URL=""
TRAIT=""
THREAD=5
NEALE_MANIFEST=~/resource/nealelab/ukb31063_ldsc_sumstat_manifest_aws_sep2022.tsv
GeneATLAS_MANIFEST=inputs/Traits_Table_GeneATLAS.csv


while [[ $# -gt 0 ]]; do
  case "$1" in
    --url)       URL="$2";       shift 2;;
    --trait)     TRAIT="$2";     shift 2;;
    --thread)    THREAD="$2";    shift 2;;
    --out_dir)   OUT_DIR="$2";   shift 2;;
    *) echo "Unknown option: $1";  exit 1;;
  esac
done


# neale lab ukbb
if [[ "$URL" == *"1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU"* ]]; then
    scripts/dl_nealeUKBB_from_a_trait.sh "$TRAIT" $NEALE_MANIFEST $OUT_DIR \
      | awk -F"\t" -v url="$URL" 'BEGIN{OFS="\t"} {print $1, url, $2, substr($0, index($0,$3))}'
    exit 0
fi


# openGWAS
if [[ "$URL" == *"gwas.mrcieu.ac.uk"* ]]; then
    source ~/.bash_profile
    scripts/dl_openGWAS.sh $URL $OUT_DIR $THREAD \
      | awk -F"\t" -v trait="$TRAIT" 'BEGIN{OFS="\t"} {print $1, $2, trait, substr($0, index($0,$3))}'
    exit 0
fi


# GeneATLAS
if [[ "$URL" == *"geneatlas.roslin.ed.ac.uk"* ]]; then
    scripts/dl_GeneATLAS_each_chr.sh "$TRAIT" $GeneATLAS_MANIFEST $OUT_DIR \
      | awk -F"\t" -v url="$URL" 'BEGIN{OFS="\t"} {print $1, url, $2, substr($0, index($0,$3))}'
    exit 0
fi


# fastGWAを含むとき
if [[ "$URL" == *"fastGWA"* ]]; then
    scripts/dl_fastGWA.sh $URL $OUT_DIR $THREAD \
      | awk -F"\t" -v trait="$TRAIT" 'BEGIN{OFS="\t"} {print $1, $2, trait, substr($0, index($0,$3))}'
    exit 0
fi


# .gzや.bgzで終わるとき
if [[ "$URL" =~ \.(gz|bgz)$ ]]; then
    f=$(basename $URL)
    if aria2c --auto-file-renaming=false --continue=true -q -s $THREAD -x $THREAD -j 1 -d $OUT_DIR $URL; then
        echo -e "0\t$URL\t$TRAIT\t$f"
        exit 0
    else
        echo -e "1\t$URL"
        exit 1
    fi
fi


# gwas catalog
if [[ "$URL" == *"GCST"* ]]; then
    scripts/dl_gwasCatalog.sh $URL $OUT_DIR $THREAD \
      | awk -F"\t" -v trait="$TRAIT" 'BEGIN{OFS="\t"} {print $1, $2, trait, substr($0, index($0,$3))}'
    exit 0
fi


# 上記のパターンに該当せず
echo -e "1\t$URL\t$TRAIT"
exit 1
