#!/usr/bin/env python3
import argparse
import pandas as pd
import tempfile
import gwaslab as gl
import sys


def main():
    parser = argparse.ArgumentParser(
            description="qtl nominalファイルをphenotypeでフィルターして、rsIDをつける、hg38のみに対応しているので注意")
    parser.add_argument("--qtl_nominal", type=str, required=True)
    parser.add_argument("--phen", type=str, required=True, help = 'ENSGxxx')
    parser.add_argument("--rsid_ref_vcf", type=str, default="/home/khamanaka/.gwaslab/GCF_000001405.40.gz")
    parser.add_argument("--rsid", type=str, required=True)
    parser.add_argument("--window", type=int, default=1000, help="kb")
    parser.add_argument("--threads", type=int, default=5)
    parser.add_argument("--out", type=str, required=True)
    args = parser.parse_args()

    ## ファイルダウンロードを初回はする
    #gl.download_ref("1kg_dbsnp151_hg38_auto")

    # parquet を読み込む
    # CHR, POS, etcを用意しないとassign_rsidが動かない、、、
    df = pd.read_parquet(args.qtl_nominal).\
        query('phenotype_id == @args.phen')

    split = df["variant_id"].str.split(":", expand=True)
    df["CHR"] = split[0].str.replace("^chr", "", case=False, regex=True)
    df["POS"] = split[1].astype(int)
    df["NEA"] = split[2]
    df["EA"]  = split[3]

    # 一時ファイルを作成（削除される）
    with tempfile.NamedTemporaryFile(suffix=".txt", mode="w", delete=False) as tmp:
        tmp_path = tmp.name
        df.to_csv(tmp_path, sep="\t", index=False)

    ss = gl.Sumstats(
        tmp_path,
        snpid="variant_id",
        chrom="CHR",
        pos="POS",
        ea="EA",
        nea="NEA",
        eaf="af",
        beta="slope",
        se="slope_se",
        p="pval_nominal"
        #n="N",
    )
    ss.basic_check()

    ss.assign_rsid(
        ref_rsid_tsv = gl.get_path("1kg_dbsnp151_hg38_auto"),
        ref_rsid_vcf = args.rsid_ref_vcf,
        chr_dict = gl.get_number_to_NC(build="38"),
        n_cores = args.threads)

    ss.filter_flanking_by_id(args.rsid, windowsizekb=args.window, inplace=True)

    ss.data.\
        rename(columns={'rsID': 'rsid', 'P': 'pval'})\
        [['rsid', 'pval']].\
        dropna().\
        to_csv(args.out, index=False, sep='\t')


if __name__ == "__main__":
    main()
