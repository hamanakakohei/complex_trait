#!/usr/bin/env python3
import argparse
import sys
import re
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract sample size info from GWAS Catalog studies table"
    )
    p.add_argument(
        "--gwa-sumsta",
        required=True,
        help="GWAS summary stats filename"
    )
    p.add_argument(
        "--traits-table",
        required=True,
        help="gwas_catalog-v1.0.3.1-studies-*.tsv"
    )
    return p.parse_args()


def extract_gcst_id(fname: str) -> str:
    m = re.search(r"(GCST\d+)", fname)
    if not m:
        sys.exit(f"GCST ID not found in GWA_SUMSTA: {fname}")
    return m.group(1)


def main():
    args = parse_args()

    gcst_id = extract_gcst_id(args.gwa_sumsta)

    try:
        df = pd.read_csv(
            args.traits_table,
            sep="\t",
            dtype=str
        )
    except Exception as e:
        sys.exit(f"Failed to read catalog file: {e}")

    required_cols = {
        "STUDY ACCESSION",
        "INITIAL SAMPLE SIZE",
        "REPLICATION SAMPLE SIZE"
    }
    if not required_cols.issubset(df.columns):
        sys.exit(f"Catalog file must contain columns: {required_cols}")

    hit = df.loc[df["STUDY ACCESSION"] == gcst_id]

    if hit.empty:
        sys.exit(f"GCST ID '{gcst_id}' not found in catalog")

    if len(hit) > 1:
        sys.exit(f"Multiple entries found for GCST ID '{gcst_id}'")

    row = hit.iloc[0]

    init_n = row["INITIAL SAMPLE SIZE"]
    repl_n = row["REPLICATION SAMPLE SIZE"]

    # join with ;
    if pd.isna(repl_n) or repl_n == "":
        N = init_n
    else:
        N = f"{init_n};{repl_n}"

    print('TRAIT_TYPE="check"')
    print(f'N="{N}"')
    print('N_CASE="check"')
    print('N_CONTROL="check"')


if __name__ == "__main__":
    main()
