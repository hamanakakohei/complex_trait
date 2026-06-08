#!/usr/bin/env python3
import argparse
import sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract trait metadata from GeneATLAS Traits table"
    )
    p.add_argument(
        "--gwa-sumsta",
        required=True,
        help="e.g. imputed.allWhites.1727-0.0.chr"
    )
    p.add_argument(
        "--traits-table",
        required=True,
        help="Traits_Table_GeneATLAS.csv"
    )
    p.add_argument(
        "--default-n",
        type=int,
        default=452264,
        help="Default sample size for quantitative traits (default: 452264)"
    )
    return p.parse_args()


def extract_key(gwa: str) -> str:
    if not gwa.startswith("imputed.allWhites.") or not gwa.endswith(".chr"):
        sys.exit(f"Unexpected GWA_SUMSTA format: {gwa}")
    return gwa.removeprefix("imputed.allWhites.").removesuffix(".chr")


def main():
    args = parse_args()

    key = extract_key(args.gwa_sumsta)

    # CSV 読み込み（"" やカンマ対策）
    try:
        df = pd.read_csv(args.traits_table)
    except Exception as e:
        sys.exit(f"Failed to read CSV: {e}")

    if "key" not in df.columns or "# Cases" not in df.columns:
        sys.exit("CSV must contain 'key' and '# Cases' columns")

    hit = df.loc[df["key"] == key]

    if hit.empty:
        sys.exit(f"Key '{key}' not found in traits table")

    if len(hit) > 1:
        sys.exit(f"Multiple entries found for key '{key}'")

    cases = str(hit.iloc[0]["# Cases"]).strip()

    if cases == "-":
        TRAIT_TYPE = "quant"
        N = args.default_n
        N_CASE = "NA"
        N_CONTROL = "NA"
    else:
        TRAIT_TYPE = "cc"
        N = "check"
        N_CASE = cases
        N_CONTROL = "check"

    # bash eval 用
    print(f'TRAIT_TYPE="{TRAIT_TYPE}"')
    print(f'N="{N}"')
    print(f'N_CASE="{N_CASE}"')
    print(f'N_CONTROL="{N_CONTROL}"')


if __name__ == "__main__":
    main()
