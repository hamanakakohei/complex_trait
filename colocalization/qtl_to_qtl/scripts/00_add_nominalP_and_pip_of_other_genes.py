#!/usr/bin/env python3
#
# newGene eqtlが他の遺伝子との
# - nominal Pがどうか？（all signifにあるか？）
# - susie cs95にあるか？全csでのmaxPIPがどうか？
# を出すのだが、この他の遺伝子には別のnewGeneが含まれるかもしれないことに注意
#
# to do:
# - "new gene"という言い方でなくて、ターゲット遺伝子か着目遺伝子みたいな一般的な名前に変えたい
import argparse
import pandas as pd


def drop_overlap(df_left, df_right, cols):
    """
    df_left から df_right に cols がすべて一致する行を除外する。
    """
    merged = df_left.merge(
        df_right[cols].drop_duplicates(),
        on=cols,
        how="left",
        indicator=True,
    )
    return merged[merged["_merge"] == "left_only"].drop(columns=["_merge"])


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--new_gene_qtl", required=True)
    parser.add_argument("--tensorqtl_all_signif_res", required=True)
    parser.add_argument("--tensorqtl_susie_res", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    newGene_eqtl = pd.read_csv(args.new_gene_qtl, sep="\t", dtype={"variant_id": str})

    # newGene_eqtl と完全一致する (phenotype_id, variant_id) の行は除外
    all_signif = pd.read_csv(
        args.tensorqtl_all_signif_res,
        sep="\t",
        usecols=["phenotype_id", "variant_id", "pval_nominal"]
      )

    all_signif = drop_overlap(
        df_left = all_signif,
        df_right = newGene_eqtl,
        cols = ["phenotype_id", "variant_id"]
      ).\
      rename(
        columns={
          "phenotype_id": "other_gene",
          "pval_nominal": "other_gene_P"
        }
      )

    # phenotype_id × variant_id ごとの最大 pip を取得
    # newGene_eqtl と重複する (phenotype_id, variant_id) の cs95 行は除外
    cs95 = pd.read_parquet(args.tensorqtl_susie_res)
    cs95_maxpip = (
        cs95.groupby(["phenotype_id", "variant_id"], as_index=False)["pip"]
        .max()
        .rename(columns={"pip": "max_pip"})
    )

    cs95_maxpip = drop_overlap(
        df_left=cs95_maxpip,
        df_right=newGene_eqtl,
        cols = ["phenotype_id", "variant_id"]
      ).\
      rename( columns={"phenotype_id": "other_gene"} )

    merged = newGene_eqtl.merge(all_signif, on="variant_id", how="left")
    merged = merged.merge(cs95_maxpip, on=["variant_id", "other_gene"], how="left")

    merged.to_csv(args.out, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
