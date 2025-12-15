#!/usr/bin/env python3
import pandas as pd
import sys
from pathlib import Path
import argparse


def parse_args():
    p = argparse.ArgumentParser(description="Bin-wise variant counting by CausalDB PIP.")
    p.add_argument("--qtl", type=Path, required=True, 
                   help="FDR1% permute結果をrsIDにしたもの（遺伝子クラスやsnpかどうかで事前にフィルターしたせいでわかりにくくなったかも、、、")
    p.add_argument("--causaldb_pip", type=Path, required=True, help="バリアントIDとPIPの2列のファイル")
    p.add_argument("--bin_table", type=Path, required=True, help="bin_name、bin_min、bin_maxの3列からならファイル")
    p.add_argument("--qtl_analysis_variants", type=Path, required=True,
                   help="qtl解析で使ったバリアントリストを与えると、それとcausaldbでかぶっているバリアントのみ考慮する")
    p.add_argument("--prefix", type=str, required=True, help="Output prefix.")
    return p.parse_args()


def main():
    args = parse_args()

    qtl = pd.read_table(args.qtl, sep="\t", dtype=str)
    causaldb = pd.read_table(args.causaldb_pip, sep="\t")
    bins = pd.read_table(args.bin_table, sep="\t")

    # バリアントリストに重複があるのはおかしいので確認する
    if qtl["rsid"].duplicated().any():
        sys.exit("[ERROR] Duplicate rsid values found in QTL file.")
    if causaldb["rsid"].duplicated().any():
        sys.exit("[ERROR] Duplicate rsid values found in CausalDB PIP file.")

    # "全バリアント"（"causaldb"）をqtl解析対象のバリアントにしぼる
    qtl_all_vars = pd.read_table(args.qtl_analysis_variants, sep="\t")[['rsid']]
    causaldb = pd.merge(causaldb, qtl_all_vars, on="rsid", how="inner")

    # QTLで、causaldb内にあるもの（pipでフィルターしない、0.0とかもある）リストを出力しておく
    df = pd.merge(qtl, causaldb, on="rsid", how="inner")
    out_variants_path = args.prefix + ".txt.gz"
    df.to_csv(out_variants_path, sep="\t", index=False)

    # === Ensure numeric PIP ===
    try:
        df["pip"] = pd.to_numeric(df["pip"], errors="raise")
    except Exception:
        sys.exit("[ERROR] PIP column contains non-numeric values.")

    # === Count variants in each bin ===
    # "causaldb"が全バリアントで、"df"がQTLかぶり
    results = []
    for _, row in bins.iterrows():
        name, minv, maxv = row["bin_name"], float(row["bin_min"]), float(row["bin_max"])
        total =    ((causaldb["pip"] >= minv) & (causaldb["pip"] <= maxv)).sum()
        n_in_bin = ((df["pip"]       >= minv) & (df["pip"]       <= maxv)).sum()
        n_out_bin = total - n_in_bin
        results.append((name, int(n_in_bin), int(n_out_bin)))

    # === Save results ===
    out_df = pd.DataFrame(results, columns=["bin", "n_in_bin", "n_out_bin"])
    out_count_path = args.prefix + ".bin_count.txt"
    out_df.to_csv(out_count_path, sep="\t", index=False)


if __name__ == "__main__":
    main()

