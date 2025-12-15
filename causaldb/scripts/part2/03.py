#!/usr/bin/env python3
import argparse
import pandas as pd


def extract_chr(variant_id):
    """
    例： chr20:1000:G:C → 20
    """
    try:
        chrom = variant_id.split(":")[0]
        return chrom.replace("chr", "")
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(description="locuscompare用の入力ファイルを作る")
    parser.add_argument("--in1", required=True)
    parser.add_argument("--in2", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    A = pd.read_csv(args.in1, sep="\t", dtype=str)
    B = pd.read_csv(args.in2, sep="\t", dtype=str)

    merged = A.merge(B, on=["url", "trait"], how="left")
    merged["chr"] = merged["variant_id"].apply(extract_chr)

    # 列順が大切、locuscompareスクリプト用に
    desired_cols = ["rsid", "chr", "phenotype_id", "f", "fmt", "popu", "url", "trait"]
    merged = merged[desired_cols]
    merged.drop_duplicates().to_csv(args.out, sep="\t", index=False)

    # ---- NA 行の抽出 ----
    na_mask = merged.isna().any(axis=1)
    merged = merged[na_mask]
    merged.to_csv("tmp.na.txt", sep="\t", index=False)


if __name__ == "__main__":
    main()
