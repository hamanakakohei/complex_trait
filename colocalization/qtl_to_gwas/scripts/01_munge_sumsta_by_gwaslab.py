#!/usr/bin/env python3
#
# この設計思想は難しくて、args.fmtで処理を分岐させた方がシンプルな気もするが、
# 同じfmt内でも複数パターンあるから、fmtで分岐させずに列の存在で分岐させた方が汎用的な気もするが、
# それだとはっきりしなくて頭が疲れる、、、
# なので、args.fmtで分岐させることにする
#
# liftoverしてからrsid付けるかその逆かだが、rsidを付ける前にSNPIDを付けることになるか？であればliftoverしてからでないと、、
# 複雑で良ければ、hg19のままSNPIDつける--> rsidつける --> SNPIDけす --> liftover --> SNPIDつける、というのもできるが、、
import argparse
import gwaslab as gl
import pandas as pd
import gzip
import sys
import math


def main():
    parser = argparse.ArgumentParser(description="指定したrsidを中心として、指定したwindowサイズに、sumstaを限定して保存する")
    parser.add_argument("--sumsta", type=str, required=True)
    parser.add_argument("--fmt", required=True,
                        help="sumstaフォーマット、参考；https://cloufield.github.io/gwaslab/tutorial_3.6/；https://github.com/Cloufield/formatbook")
    parser.add_argument("--rsid", type=str, default=None, help="指定したrsidがsumsta内になければエラーにする")
    parser.add_argument("--threads", type=int, default=5)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # infer_buildするとSTATUSがNaNになることがあり、一度inferしたらもう一度読み直す
    ss = gl.Sumstats(args.sumsta, fmt=args.fmt, verbose=True)
    cols = ss.data.columns

    if "CHR" not in cols or "POS" not in cols:
      if 'SNPID' in cols:
        ss.fix_id(fixchrpos=True, fixeanea =True, fixprefix=True)
      else:
        # rsIDのみなら勝手にhg38にする
        ss.rsid_to_chrpos( path = gl.get_path("1kg_dbsnp151_hg38_auto"))

    ss.infer_build()
    build = ss.meta['gwaslab']['genome_build']

    ss = gl.Sumstats(args.sumsta, fmt=args.fmt, verbose=True)
    ss.data = ss.data.dropna(axis=1, how='all') # 全行がNAの列を落とす
    cols = ss.data.columns

    print(1)
    print(ss.data)
    print(ss.data.columns)

    if args.fmt == "NealeLab":
      ss.data["EAF"] = ss.data["AC"] / (ss.data["N"] * 2)
      ss.fix_id(fixchrpos=True, fixeanea =True, fixprefix=True)

    if "CHR" not in cols or "POS" not in cols:
      if 'SNPID' in cols:
        ss.fix_id(fixchrpos=True, fixeanea =True, fixprefix=True)
      else:
        # rsIDのみなら勝手にhg38にする
        ss.rsid_to_chrpos( path = gl.get_path("1kg_dbsnp151_hg38_auto"))

    cols = ss.data.columns

    print(2)
    print(ss.data)
    print(ss.data.columns)

    # これgwaslabのバグな気がする、こうしないとPが作れない？
    if "MLOG10P" in cols:
      ss.data["MLOG10P"] = pd.to_numeric(ss.data["MLOG10P"], errors="coerce")

    if "P" not in cols:
      ss.fill_data(to_fill=["P"])

    if "BETA" not in cols and "OR" in cols:
      ss.data["BETA"] = ss.data["OR"].map(lambda x: math.log(x))

    if "VAR_BETA" not in cols and "SE" in cols:
      ss.data["SE"] = pd.to_numeric(ss.data["SE"], errors="coerce")
      ss.data["VAR_BETA"] = ss.data["SE"] ** 2

    # NAがあるとCHR列がfloatとして読まれて、.basic_check()で捨てられるので
    ss.data = ss.data.dropna(subset='CHR')
    if pd.api.types.is_float_dtype(ss.data["CHR"]):
      ss.data["CHR"] = ss.data["CHR"].astype(str)

    # OR列などにNAがあると、.basic_check()で捨てられるので
    excl_cols = ['OR','OR_95L','OR_95U']
    excl_cols = [c for c in excl_cols if c in ss.data.columns]
    ss.data = ss.data.drop(columns=excl_cols)

    # ここでP列が1e-300以下なら1e-300にするとかしたほうがいいかも

    # checkは列が揃ってからの方が良いのでここで
    ss.basic_check()

    print(5)
    print(ss.data)
    print(ss.data.columns)

    #if args.fmt == "gwascatalog_hm" or "fastgwa" or "ssf":
    # POSをhg38にする（INDELはもはやあてにならない？
    if build == "19":
      ss.liftover(n_cores=args.threads, from_build="19", to_build="38")

      # NEAがrefか確認して、必要ならstatsのプラマイを変える？？
      ss.check_ref(gl.get_path("ucsc_genome_hg38"))
      ss.flip_allele_stats()

      if "SNPID" in cols:
        ss.data = ss.data.drop(columns='SNPID')

      ss.fix_id(fixid=True, forcefixid=True)

    ss.check_ref(gl.get_path("ucsc_genome_hg38"))
    ss.flip_allele_stats()

    if "SNPID" not in cols:
      ss.fix_id(fixid=True, forcefixid=True)

    if "rsID" not in cols:
      if 'SNPID' not in cols:
        ss.fix_id(fixid=True, forcefixid=True)

      # SNPID列にrsIDが混ざっているならそれをrsID列としても使う
      if ss.data['SNPID'].fillna('').str.startswith('rs').any():
        ss.data['rsID'] = ss.data['SNPID']
      else:
        ss.fix_id(fixid=True, forcefixid=True)
        ss.assign_rsid(
          ref_rsid_tsv = gl.get_path("1kg_dbsnp151_hg38_auto"),
          chr_dict = gl.get_number_to_NC(build = "38"),
          n_cores = args.threads
        )

    # 興味のあるrsIDがなければ、ポジションから付け直す
    #if args.rsid is not None and args.rsid not in ss.data['rsID']:
    if args.rsid and not ss.data["rsID"].eq(args.rsid).any():
      ss.assign_rsid(
        ref_rsid_vcf = gl.get_path("dbsnp_v156_hg" + build),
        chr_dict = gl.get_number_to_NC(build = build),
        n_cores = args.threads
      )

      # それでもなければエラー
      if not ss.data["rsID"].eq(args.rsid).any():
        raise ValueError(f"args.rsid not found in the sumsta")

    print(6)
    print(ss.data)
    print(ss.data.columns)

    # SNPIDを作りなおす、本当は必要ないと思うがSNPIDを間違ってつけるssfを使ってしまっているから、、、
    ss.data = ss.data.drop(columns='SNPID')
    ss.fix_id(fixid=True, forcefixid=True)
    #ss.fix_id(fixsep=True)
    ss.data.to_csv(args.out, index=False, sep='\t')


if __name__ == "__main__":
    main()
