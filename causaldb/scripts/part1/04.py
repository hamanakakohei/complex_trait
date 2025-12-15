#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Filter CausalDB fine-mapping results and merge with metadata.")
    parser.add_argument("--causaldb_meta", required=True, type=Path, help="meta.txt")
    parser.add_argument("--causaldb_cs", required=True, type=Path, help="01.pyでフィルター後やつか元のcredible_set.txt.gz")
    parser.add_argument("--causaldb_fm_method", required=True,
                        choices=["pip", "abf", "finemap", "paintor", "caviarbf", "susie", "polyfun_finemap", "polyfun_susie"],
                        help="上で01.pyでフィルター後のファイルを指定した時は「pip」、credible.txt.gzの時はFM手法列名")
    parser.add_argument("--causaldb_pip_thr", required=True, type=float)
    parser.add_argument("--target_rsid_file", required=True, type=Path, help="02.pyの結果、rsid列を使いたい")
    parser.add_argument("--qtl_meta", required=True, type=Path, help="tensorqtl indepファイルでrsidにqtl情報を加える")
    parser.add_argument("--rsid_varid_map", required=True, help="rsidとvariant_idの2列")
    parser.add_argument("--out1", required=True, type=Path)
    parser.add_argument("--out2", required=True, type=Path)
    args = parser.parse_args()

    df_cs = pd.read_csv(args.causaldb_cs, sep="\t")
    df_causaldb_meta = pd.read_csv(args.causaldb_meta, sep="\t")
    df_target = pd.read_csv(args.target_rsid_file, sep="\t")
    df_qtl_meta = pd.read_csv(args.qtl_meta, sep="\t")[['phenotype_id', 'variant_id', 'rank']]
    df_rsid_varid = pd.read_csv(args.rsid_varid_map, sep="\t")[["rsid", "variant_id"]]


    # causaldbをターゲットrsidとpip値でフィルターする
    rsids_oi = set(df_target["rsid"].astype(str))

    df_filtered = (
        df_cs[
            (df_cs["rsid"].astype(str).isin(rsids_oi)) &
            (df_cs[args.causaldb_fm_method] >= args.causaldb_pip_thr)
        ][["rsid", "meta_id", args.causaldb_fm_method]]
        .drop_duplicates()
    )


    # causaldbメタ情報を付ける
    df_out = df_filtered.merge(df_causaldb_meta, on="meta_id", how="left")


    # qtlメタ情報を付ける
    df_qtl_with_rsid = df_qtl_meta.merge(
        df_rsid_varid,
        on="variant_id",
        how="inner"
    )

    df_out = df_out.merge(
        df_qtl_with_rsid,
        on="rsid",
        how="left"
    )


    # NA埋めて保存
    df_out = df_out.sort_values(by=["rsid", args.causaldb_fm_method])
    args.out1.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.out1, sep="\t", index=False, na_rep="NA")
    df_out[['url', 'trait']].drop_duplicates().to_csv(args.out2, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
