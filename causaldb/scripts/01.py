#!/usr/bin/env python3
from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Filter and summarize CausalDB credible set table.")
    p.add_argument("--causaldb", type=Path, required=True)
    p.add_argument("--method", type=str, choices=["abf", "finemap", "paintor", "caviarbf", "susie", "polyfun_finemap", "polyfun_susie"])
    p.add_argument("--snp_only", action="store_true")
    p.add_argument("--cs95_only", action="store_true", help="causaldb側のCSのPIP合計が95以上にならないときがあるので")
    p.add_argument("--primary_only", action="store_true", help="primary signalかどうか")
    p.add_argument("--p_thr", type=float, default=5e-8, help="causaldb側のフィルター")
    p.add_argument("--prefix", type=str, required=True)
    return p.parse_args()

def main():
    args = parse_args()

    df = pd.read_csv(args.causaldb, sep="\t", compression="infer", low_memory=False)
    df = df.rename(columns={args.method: "pip"})
    df["bp"] = df["bp"].astype(str)
    df["chr"] = df["chr"].astype(str)
    df["variant_id"] = df["chr"] + "_" + df["bp"] + "_" + df["nea"] + "_" + df["ea"]

    # pip列が-1やNaNは0にする、かつ数値型を確認する
    try:
        df["pip"] = pd.to_numeric(df["pip"], errors="raise")
    except Exception:
        sys.exit("[ERROR] PIP column contains non-numeric values.")

    df["pip"] = df["pip"].replace(-1, 0).fillna(0)

    # ea, nea列にA,C,G,T以外のもじがあるとエラー
    for allele_col in ["ea", "nea"]:
        invalid_mask = ~df[allele_col].astype(str).str.fullmatch(r"[ACGT]+", na=False)
        if invalid_mask.any():
            sys.exit(f"[ERROR] Invalid characters found in '{allele_col}' column.")

    if args.cs95_only:
        grp = df.groupby(["meta_id", "block_id", "primary"], observed=True)["pip"].sum()
        good = grp[grp > 0.95].reset_index()[["meta_id", "block_id", "primary"]]
        df = df.merge(good, on=["meta_id", "block_id", "primary"], how="inner")

    if args.snp_only:
        df = df[(df["ea"].str.len() == 1) & (df["nea"].str.len() == 1)]

    if args.primary_only:
        df = df[df["primary"].astype(str) == "1"]

    # pip列が-1やNaNは0にする、かつ数値型を確認する
    try:
        df["p"] = pd.to_numeric(df["p"], errors="raise")
    except Exception:
        sys.exit("[ERROR] Column 'p' contains non-numeric values.")

    df = df[df["p"] < args.p_thr]

    # そのまま出力するのと、バリアントごとにmax PIPを出すのと
    out_main = args.prefix + ".txt.gz"
    df.to_csv(out_main, sep="\t", index=False)

    out_rsid = args.prefix + ".rsid_156.maxPIP.txt.gz"
    df.groupby("rsid", observed=True)["pip"].max().reset_index().to_csv(out_rsid, sep="\t", index=False)

    out_var = args.prefix + ".varid_hg19.maxPIP.txt.gz"
    df.groupby("variant_id", observed=True)["pip"].max().reset_index().to_csv(out_var, sep="\t", index=False)


if __name__ == "__main__":
    main()
