#!/usr/bin/env python3
import sys
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract trait metadata from Pan-UK Biobank phenotype manifest"
    )
    parser.add_argument(
        "--gwa-sumsta",
        required=True,
        help="GWA summary statistics filename (must match 'filename' column)",
    )
    parser.add_argument(
        "--traits-table",
        required=True,
        help="Path to Pan-UK Biobank phenotype manifest Excel file or phenotype_manifest sheet .txt",
    )
    parser.add_argument(
        "--sheet",
        default="phenotype_manifest",
        help="Sheet name in the Excel file (default: phenotype_manifest)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    gwa = args.gwa_sumsta

    # Excel 読み込み
    try:
        if args.traits_table.endswith((".xlsx", ".xls")):
            df = pd.read_excel(args.traits_table, sheet_name=args.sheet)
        else:
            df = pd.read_table(args.traits_table)
    except Exception as e:
        sys.exit(f"Failed to read Excel file: {e}")

    # 必須列チェック
    required_cols = {
        "filename",
        "trait_type",
        "n_cases_EUR",
        "n_controls_EUR",
    }
    missing = required_cols - set(df.columns)
    if missing:
        sys.exit(f"Missing columns in Excel: {missing}")

    # filename 完全一致
    hit = df.loc[df["filename"] == gwa]

    if hit.empty:
        sys.exit(f"{gwa} not found")

    if len(hit) > 1:
        sys.exit(f"Multiple entries found for {gwa}")

    row = hit.iloc[0]

    TRAIT_TYPE = row["trait_type"]
    N_CASE = row["n_cases_EUR"]
    N_CONTROL = row["n_controls_EUR"]

    if TRAIT_TYPE in {"categorical", "icd10", "phecode", "prescriptions"}:
        TRAIT_TYPE = "cc"
    elif TRAIT_TYPE in {"biomarkers", "continuous"}:
        TRAIT_TYPE = "quant"
    else:
        sys.exit(f"Unknown variable_type: {TRAIT_TYPE}")


    # NA 対策
    if pd.isna(N_CASE):
        N_CASE = "NA"
    else:
        N_CASE = int(N_CASE)

    if pd.isna(N_CONTROL):
        N_CONTROL = "NA"
    else:
        N_CONTROL = int(N_CONTROL)

    if N_CASE == "NA" and N_CONTROL == "NA":
        N = "NA"
    elif N_CASE == "NA":
        N = N_CONTROL
    elif N_CONTROL == "NA":
        N = N_CASE
    else:
        N = N_CASE + N_CONTROL

    # bash eval 用出力
    print(f'TRAIT_TYPE="{TRAIT_TYPE}"')
    print(f'N="{N}"')
    print(f'N_CASE="{N_CASE}"')
    print(f'N_CONTROL="{N_CONTROL}"')


if __name__ == "__main__":
    main()
