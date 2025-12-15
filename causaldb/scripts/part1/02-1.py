#!/usr/bin/env python3
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="QTL indepをフィルターしつつrsidにする")
    parser.add_argument("--qtl_indep", required=True, help="tensorqtl indepの結果")
    parser.add_argument("--gene_class_map", required=True, help="gene_idとstructural_categoryの2列")
    parser.add_argument("--gene_class", required=True, nargs='+',
                        help="Gene class filter (space-separated list). "
                             "Special keywords: 'new', 'known'")
    parser.add_argument("--rsid_varid_map", required=True, help="rsidとvariant_idの2列")
    parser.add_argument("--primary_only", action='store_true', help="tensorqtl indepのrank列==1")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    df_qtl = pd.read_csv(args.qtl_indep, sep="\t")
    

    # gene_class フィルタリング
    new_set = {"intergenic", "genic_intron", "antisense"}
    selection = set(args.gene_class)

    df_gene = pd.read_csv(args.gene_class_map, sep="\t")
    df_merged = df_qtl.merge(df_gene[['gene_id', 'structural_category']], left_on='phenotype_id', right_on='gene_id', how='inner')

    if selection == {"new"}:
        df_merged = df_merged[ df_merged["structural_category"].isin(new_set)]
    elif selection == {"known"}:
        df_merged = df_merged[~df_merged["structural_category"].isin(new_set)]
    else:
        df_merged = df_merged[ df_merged["structural_category"].isin(selection)]


    # rsid mapping
    df_rsid = pd.read_csv(args.rsid_varid_map, sep="\t")
    df_merged = df_merged.merge(df_rsid[['variant_id', 'rsid']], on='variant_id', how='inner')


    # primary_only
    if args.primary_only:
        df_merged = df_merged[df_merged['rank'] == 1]
    
    df_merged[['rsid']].drop_duplicates().to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
