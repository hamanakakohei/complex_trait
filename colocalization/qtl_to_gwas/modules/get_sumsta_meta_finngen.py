#!/usr/bin/env python3
import argparse
import sys
import re
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract trait metadata from FinnGen endpoints.tsv"
    )
    p.add_argument(
        "--gwa-sumsta",
        required=True,
        help="GWAS summary stats filename (e.g. finngen_R9_L12_DERMATITISECZEMA.gz)"
    )
    p.add_argument(
        "--traits-table",
        required=True,
        help="endpoints.tsv"
    )
    return p.parse_args()


def extract_phenocode(gwa_sumsta: str) -> str:
    # remove .gz
    x = re.sub(r"\.gz$", "", gwa_sumsta)
    # remove finngen_R[digits]_ prefix
    x = re.sub(r"^finngen_R\d+_", "", x)
    return x


def main():
    args = parse_args()

    phenocode = extract_phenocode(args.gwa_sumsta)

    try:
        df = pd.read_csv(
            args.traits_table,
            sep="\t",
            dtype=str
        )
    except Exception as e:
        sys.exit(f"Failed to read endpoints.tsv: {e}")

    required_cols = {
        "phenocode", "num_cases", "num_controls"
    }
    if not required_cols.issubset(df.columns):
        sys.exit(f"endpoints.tsv must contain columns: {required_cols}")

    hit = df.loc[df["phenocode"] == phenocode]

    if hit.empty:
        sys.exit(f"phenocode '{phenocode}' not found in endpoints.tsv")

    if len(hit) > 1:
        sys.exit(f"Multiple entries found for phenocode '{phenocode}'")

    row = hit.iloc[0]

    num_cases = int(row["num_cases"])
    num_controls = int(row["num_controls"])

    N = num_cases + num_controls
    N_CASE = num_cases
    N_CONTROL = num_controls

    if num_controls == 0:
        TRAIT_TYPE = "quant"
    else:
        TRAIT_TYPE = "cc"

    # bash eval 用
    print(f'TRAIT_TYPE="{TRAIT_TYPE}"')
    print(f'N="{N}"')
    print(f'N_CASE="{N_CASE}"')
    print(f'N_CONTROL="{N_CONTROL}"')


if __name__ == "__main__":
    main()
