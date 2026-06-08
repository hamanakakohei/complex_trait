#!/usr/bin/env python3
import argparse
import sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract trait metadata from OpenGWAS CSV"
    )
    p.add_argument(
        "--gwa-sumsta",
        required=True,
        help="e.g. met-d-XL_HDL_C_pct.vcf.gz"
    )
    p.add_argument(
        "--traits-table",
        required=True,
        help="opengwas.io_*.csv"
    )
    return p.parse_args()


def extract_opengwas_id(gwa: str) -> str:
    if not gwa.endswith(".vcf.gz"):
        sys.exit(f"Unexpected GWA_SUMSTA format: {gwa}")
    return gwa.removesuffix(".vcf.gz")


def main():
    args = parse_args()

    opengwas_id = extract_opengwas_id(args.gwa_sumsta)

    try:
        df = pd.read_csv(args.traits_table)
    except Exception as e:
        sys.exit(f"Failed to read CSV: {e}")

    required_cols = {"OpenGWAS ID", "Category", "Sample Size"}
    if not required_cols.issubset(df.columns):
        sys.exit(f"CSV must contain columns: {required_cols}")

    hit = df.loc[df["OpenGWAS ID"] == opengwas_id]

    if hit.empty:
        sys.exit(f"OpenGWAS ID '{opengwas_id}' not found")

    if len(hit) > 1:
        sys.exit(f"Multiple entries found for OpenGWAS ID '{opengwas_id}'")

    row = hit.iloc[0]
    category = str(row["Category"]).strip()
    sample_size = row["Sample Size"]

    # Sample Size が空のとき
    if pd.isna(sample_size) or sample_size == "":
        N = "NA"
    else:
        N = str(int(sample_size)) #if str(sample_size).isdigit() else str(sample_size)

    if category in {"Categorical Ordered", "Continuous"}:
        TRAIT_TYPE = "quant"
        N_CASE = "NA"
        N_CONTROL = "NA"

    elif category in {"Binary", "Disease"}:
        TRAIT_TYPE = "cc"
        N_CASE = "check"
        N_CONTROL = "check"

    else:
        TRAIT_TYPE = "check"
        N_CASE = "check"
        N_CONTROL = "check"

    # bash eval 用に出力
    print(f'TRAIT_TYPE="{TRAIT_TYPE}"')
    print(f'N="{N}"')
    print(f'N_CASE="{N_CASE}"')
    print(f'N_CONTROL="{N_CONTROL}"')


if __name__ == "__main__":
    main()
