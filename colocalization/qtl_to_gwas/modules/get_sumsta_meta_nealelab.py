#!/usr/bin/env python3
import sys
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract trait metadata from Pan-UK Biobank TSV (Neale lab)"
    )
    parser.add_argument(
        "--gwa-sumsta",
        required=True,
        help="GWA summary statistics filename",
    )
    parser.add_argument(
        "--traits-table",
        required=True,
        help="Path to phenotypes.both_sexes.v2.tsv.bgz",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    gwa = args.gwa_sumsta
    tsv = args.traits_table

    # phenotype 名を GWA_SUMSTA から抽出（最初の . より前）
    phenotype = gwa.split(".", 1)[0]

    # TSV 読み込み（bgzip 対応）
    try:
        df = pd.read_csv(tsv, sep="\t", compression="gzip")
    except Exception as e:
        sys.exit(f"Failed to read TSV file: {e}")

    # 必須列チェック
    required_cols = {
        "phenotype",
        "variable_type",
        "n_cases",
        "n_controls",
    }
    missing = required_cols - set(df.columns)
    if missing:
        sys.exit(f"Missing columns in TSV: {missing}")

    # phenotype 一致行を抽出
    hit = df.loc[df["phenotype"] == phenotype]

    if hit.empty:
        sys.exit(f"{phenotype} not found in {tsv}")

    if len(hit) > 1:
        sys.exit(f"Multiple entries found for phenotype: {phenotype}")

    row = hit.iloc[0]

    # variable_type → TRAIT_TYPE 変換
    vt = row["variable_type"]

    if vt == "binary":
        TRAIT_TYPE = "cc"
    elif vt in {"continuous_irnt", "continuous_raw", "ordinal"}:
        TRAIT_TYPE = "quant"
    else:
        sys.exit(f"Unknown variable_type: {vt}")

    # n_cases / n_controls
    N = row["n_non_missing"]
    N_CASE = row["n_cases"]
    N_CONTROL = row["n_controls"]

    # NA 対策
    if pd.isna(N):
        N = "NA"
    else:
        N = int(N)

    if pd.isna(N_CASE):
        N_CASE = "NA"
    else:
        N_CASE = int(N_CASE)

    if pd.isna(N_CONTROL):
        N_CONTROL = "NA"
    else:
        N_CONTROL = int(N_CONTROL)

    # bash eval 用出力
    print(f'TRAIT_TYPE="{TRAIT_TYPE}"')
    print(f'N="{N}"')
    print(f'N_CASE="{N_CASE}"')
    print(f'N_CONTROL="{N_CONTROL}"')


if __name__ == "__main__":
    main()
