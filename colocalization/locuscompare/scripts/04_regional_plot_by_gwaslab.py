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
    parser.add_argument("--plot_build", type=str, default="38", help="このbuildで（必要ならliftoverして）プロットする")
    parser.add_argument("--gtf", type=str, required=True)
    parser.add_argument("--pop", type=str, default="EUR", help="eas, eur, sas, amr, afr")
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
    ss.infer_build()
    build = ss.meta['gwaslab']['genome_build']

    if "rsID" not in cols:
      if 'SNPID' not in cols:
        ss.fix_id(fixid=True, forcefixid=True)

      if ss.data['SNPID'].fillna('').str.startswith('rs').any():
        # SNPID列にrsIDが混ざっているならそのままrsID列として使う
        ss.data['rsID'] = ss.data['SNPID']
      else:
        # 混ざっていなくてがちのSNPIDなら、ちゃんとrsID列を作る
        ss.fix_id(fixid=True, forcefixid=True)
        ss.assign_rsid(
          ref_rsid_tsv = gl.get_path("1kg_dbsnp151_hg" + build + "_auto"),
          chr_dict = gl.get_number_to_NC(build = build),
          n_cores = args.threads)
        if args.rsid not in ss.data['rsID']:
          ss.assign_rsid(
            ref_rsid_vcf = gl.get_path("dbsnp_v156_hg" + build),
            chr_dict = gl.get_number_to_NC(build = build),
            n_cores = args.threads)


    # ここで領域をせばめる
    # 事前にフィルターされていないと次のリフトオーバーで時間かかるから
    ss.filter_flanking_by_id(args.rsid, windowsizekb=args.window, inplace=True)

    # リフトオーバー
    if build != args.plot_build:
      ss.liftover(
        from_build = build,
        to_build = args.plot_build,
        remove=True,
        n_cores=args.threads)

    CHR = ss.data['CHR'].mode()[0]
    STA = ss.data['POS'].min()
    END = ss.data['POS'].max()
    #gl.dump_pickle(ss, "tmp4.pickle", overwrite=False)

    # 引数を辞書にまとめる
    kwargs = dict(
      mode = "r",
      region = (CHR, STA, END),
      anno = "rsID",
      anno_set = [args.rsid],
      gtf_path = args.gtf,
      region_ref = args.rsid,
      region_grid = True,
      region_lead_grid = True,
      region_protein_coding = False,
      build = args.plot_build,
      title = args.out
    )

    # EA も NEA も無いときだけ vcf_path を追加
    if ("EA" in ss.data.columns) and ("NEA" in ss.data.columns):
      kwargs["vcf_path"] = gl.get_path(
        f"1kg_{args.pop.lower()}_hg{args.plot_build}"
      )

    fig, log = ss.plot_mqq(**kwargs)
    fig.savefig(args.out, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
