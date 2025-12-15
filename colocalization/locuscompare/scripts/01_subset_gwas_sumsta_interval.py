#!/usr/bin/env python3
import argparse
import gwaslab as gl
import pandas as pd
import gzip


def main():
    parser = argparse.ArgumentParser(description="指定したrsidを中心として、指定したwindowサイズに、sumstaを限定して保存する")
    parser.add_argument("--sumsta", type=str, required=True)
    parser.add_argument("--fmt", required=True,
                        help="sumstaフォーマット、参考；https://cloufield.github.io/gwaslab/tutorial_3.6/；https://github.com/Cloufield/formatbook")
    parser.add_argument("--rsid", type=str, required=True)
    parser.add_argument("--window", type=int, default=1000, help="kb")
    parser.add_argument("--threads", type=int, default=5)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    ss = gl.Sumstats(args.sumsta, fmt=args.fmt, verbose=True)
    ss.data = ss.data.dropna(axis=1, how='all')
    cols = ss.data.columns

    if "CHR" not in cols or "POS" not in cols:
      if 'SNPID' in cols:
        ss.fix_id(fixchrpos=True, fixeanea =True, fixprefix=True)
      else:
        # rsIDのみなら勝手にhg38にする
        ss.rsid_to_chrpos( path = gl.get_path("1kg_dbsnp151_hg38_auto"))

    # これgwaslabのバグな気がする、こうしないとPが作れない？
    if "MLOG10P" in cols:
      ss.data["MLOG10P"] = pd.to_numeric(ss.data["MLOG10P"], errors="coerce")

    if "P" not in cols:
      ss.fill_data(to_fill=["P"])

    # ここでP列が1e-300以下なら1e-300にするとかしたほうがいいかも

    # checkは列が揃ってからの方が良いのでここで
    ss.basic_check()

    if "rsID" not in cols:
      if 'SNPID' not in cols:
        ss.fix_id(fixid=True, forcefixid=True)

      # SNPID列にrsIDが混ざっていることがあり、その時はOK
      if not ss.data['SNPID'].fillna('').str.startswith('rs').any():
        ss.fix_id(fixid=True, forcefixid=True)
        ss.infer_build()
        build = ss.meta['gwaslab']['genome_build']
        ss.assign_rsid(
          ref_rsid_tsv = gl.get_path("1kg_dbsnp151_hg" + build + "_auto"),
          chr_dict = gl.get_number_to_NC(build = build),
          n_cores = args.threads)
        if args.rsid not in ss.data['rsID']:
          ss.assign_rsid(
            ref_rsid_vcf = gl.get_path("dbsnp_v156_hg" + build),
            chr_dict = gl.get_number_to_NC(build = build),
            n_cores = args.threads)

    ss\
      .filter_flanking_by_id(args.rsid, windowsizekb=args.window)\
      .data\
      .rename(columns=( {"rsID": "rsid"} if "rsID" in ss.data.columns else {"SNPID": "rsid"} ))\
      .rename(columns={'P': 'pval'})\
      [['rsid', 'pval']]\
      .dropna()\
      .to_csv(args.out, index=False, sep='\t')


if __name__ == "__main__":
    main()
