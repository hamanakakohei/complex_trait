#!/usr/bin/env python3
import argparse
import sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract GWAS trait metadata from manual mapping table"
    )
    p.add_argument(
        "--gwa-sumsta",
        required=True,
        help="GWAS summary statistics filename"
    )
    p.add_argument(
        "--traits-table",
        required=True,
        help="gwasatlas_trait_info_manual.txt"
    )
    return p.parse_args()


def main():
    args = parse_args()

    try:
        df = pd.read_csv(
            args.traits_table,
            sep=r"\s+",
            dtype=str
        )
    except Exception as e:
        sys.exit(f"Failed to read manual table: {e}")

    required_cols = {
        "N", "N_CASE", "N_CONTROL", "TRAIT_TYPE", "GWAS_SUMSTA"
    }
    if not required_cols.issubset(df.columns):
        sys.exit(f"Manual table must contain columns: {required_cols}")

    hit = df.loc[df["GWAS_SUMSTA"] == args.gwa_sumsta]

    if hit.empty:
        sys.exit(f"GWAS_SUMSTA '{args.gwa_sumsta}' not found in manual table")

    if len(hit) > 1:
        sys.exit(f"Multiple entries found for GWAS_SUMSTA '{args.gwa_sumsta}'")

    row = hit.iloc[0]

    N = row["N"]
    N_CASE = row["N_CASE"]
    N_CONTROL = row["N_CONTROL"]
    TRAIT_TYPE = row["TRAIT_TYPE"]

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
