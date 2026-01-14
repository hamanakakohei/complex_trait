to do
- gwas catalog
  - ダウンロードをapiでする
  - urlが間違っていてもgcst idから自動で出来るようにする
- 他のdbも既製のスクリプトがないのか？？

**GeneATLAS**
- trait表：http://geneatlas.roslin.ed.ac.uk/traits-table/
  - ここでcase controlか連続形質か、case数control数がわかるが、連続形質の時は全サンプル数がわからない、、
  - 一覧表をダウンロードできる：Traits_Table_GeneATLAS.csv
- FAQ: http://geneatlas.roslin.ed.ac.uk/frequently-asked-questions/
  -  1st release: 408,455 white-British; 2nd release: 452,264 white ancestry らしいので連続形質の時はこれを全サンプル数とするか
- sumsta
  - ダウンロード用スクリプトはhttp://geneatlas.roslin.ed.ac.uk/downloads/ で適当なtraitを入れると汎用的なものが得られる
  - 列（参照：http://www.dissect.ed.ac.uk/documentation-gwas/）
    - SNP ALLELE iscores NBETA-xxx NSE-xxx PV-xxx
    - このALLELEはrefで下のsnp情報ファイルのA2と一致する、なのでaltは下のA1がそれ
- snp情報
  - ダウンロード: http://geneatlas.roslin.ed.ac.uk/downloads/ でsumstaのとなりのstatsファイルから
  - 列
    - SNP Position A1 A2 MAF HWE-P iscore
    - このMAFは本当にminor alleleの方なのでA1 freqかA2 freqかはわからない、、

**ieu openGWAS**
- trait表: https://opengwas.io/datasets/
  - ここの「Export current results to CSV」で全traitのsample sizeがわかる
  - 下の「Search:」でtraitを検索してさらにクリックすればcase数がわかる
  - hg19 (https://opengwas.io/about)

**GWASATLAS**
- 参考：https://atlas.ctglab.nl/documentation
- trait表：https://atlas.ctglab.nl/traitDB
  - ここの下の「Search:」でtraitを検索するとcase数control数やらわかる
  - システマティックにこれらを得る方法は無さそう
- sumsta
  - GRCh37
  - BETAかORかどちらか
  - SNP列のアレル表記は「alphabetically ordered」とのことで何の意味も無し、、、
  - A1はMAFが低いのを選んでいそうなので、MAF = EAFとなる

**gwas catalog**
- 注意点：
  - hm_chromはNAのことがあり（chromosome列は値があるのに）読み込むとfloatになるときがある
  - OR等の列がNAのままだと.basic_check()でそのバリアントは除かれるので、列を消しておく
  - case control数などどこにものっていない時もあり、論文を見ないといけないかも
- GCSTxxx.h.tsv.gz
  - **chromosome, base_pair_location**, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, **p_value**, variant_id, hm_coordinate_conversion, hm_code, **rsid**
    - Direction, HetChiSq, HetDf, HetISq, HetPVal, 
    - REF, n, INFO
    - markername, freqse, minfreq, maxfreq, direction, hetisq, hetchisq, hetdf, hetpval, cases, effective_cases, n, odds_ratio
    - markername, freqse, minfreq, maxfreq, direction, hetisq, hetchisq, hetdf, hetpval, cases, effective_cases, n, meta_analysis, odds_ratio
    - n, N_studies
    - GENPOS, INFO, CHISQ_LINREG, P_LINREG, CHISQ_BOLT_LMM_INF, P_BOLT_LMM_INF, CHISQ_BOLT_LMM
    - variant_id_hg19, base_pair_location_grch38
- xxx-GCSTxxx-<EFO or GO>_xxx.h.tsv.gz
  - hm_variant_id, **hm_rsid**, **hm_chrom, hm_pos**, hm_other_allele, hm_effect_allele, hm_beta, hm_odds_ratio, hm_ci_lower, hm_ci_upper, hm_effect_allele_frequency, hm_code, variant_id, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, **p_value**, odds_ratio, ci_lower, ci_upper
    - n
    - n, snp.1, info
    - n, uniqid, af, info, zval
    - range

**Neale lab**
- 参照：
  - FAQ: https://www.nealelab.is/uk-biobank/faq
  - マニフェスト: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=227859291#gid=227859291
    - 主なファイルの列の説明
      - variants.tsv.bgz
      - phenotypes.{both_sexes,female,male}.tsv.bgz
        - EAFはAC / (n_complete_samples * 2)と自分で計算するのがよさそう 
      - <phenotype_code>.gwas.imputed_v3.{both_sexes,female,male}.tsv.bgz
    - 他にも色々なメタデータファイルあり：
      - Updated v2 phenotype summary file（phenotypes.both_sexes.v2.tsv.bgz）
        - variable_type列：binary, continuous_irnt, continuous_raw, ordinal
        - n_non_missing, n_controls, n_cases列
      - List of phenotypes with updated sample counts in v2 phenotype summary file
      - List of variants used in GWAS, with annotations 
      - Summary of biomarker phenotypes
  - binary phenotypeもlinear regressionしている？：https://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas
- ordinal
  - どう解析した？？？
- biomarkers
  - https://www.nealelab.is/blog/2019/9/16/biomarkers-gwas-results
- xxx_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
  - **variant**, minor_allele, minor_AF, low_confidence_variant, n_complete_samples, AC, ytx, beta, se, tstat, **pval**

**pan-ukbb**
- 参照：
  - per-phenotype file: https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files/index.html
    - "continuous" or "biomarkers" --> quantitative
    - "prescriptions", "icd10", "phecode" or "categorical" --> binary
  - manifest: https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit?gid=1450719288#gid=1450719288
    - trait_type列：continuous, biomarkers, prescriptions, icd10, phecode, categoricalのどれか
    - n_cases_full_cohort_both_sexes（これのn_controls版が無い）, n_cases_EUR, n_controls_EUR
    - filename列
  - LD matrices: https://pan.ukbb.broadinstitute.org/docs/hail-format#exporting-a-ld-matrix-to-a-flat-file
  - https://pan-dev.ukbb.broadinstitute.org/downloads
- biomarkers-xxx-both_sexes-irnt.tsv.bgz 
  - **chr, pos**, ref, alt
    - af_<meta_hq or meta or each_of_6_pops>
    - beta_<meta_hq or meta or each_of_6_pops>
    - se_<meta_hq or meta or each_of_6_pops>
    - **neglog10_pval_**<meta_hq or meta or each_of_6_pops>
    - neglog10_pval_heterogeneity_hq, neglog10_pval_heterogeneity
    - low_confidence_<each_of_6_pops>
- categorical-xxx-both_sexes-xxx.tsv.bgz
  - **chr, pos**, ref, alt
    - af_cases_<meta or each_of_6_pops>
    - af_controls_<meta or each_of_6_pops>
    - beta_<meta or each_of_6_pops>
    - se_<meta or each_of_6_pops>
    - **neglog10_pval_**<meta or each_of_6_pops>
    - neglog10_pval_heterogeneity
    - low_confidence_<each_of_6_pops>
- continuous-xxx--both_sexes-<"" or "auto_medadj" or "combined_medadj">_<irnt or raw>.tsv.bgz
  - xxx_meta_hqはあったり無かったりなので、xxx_metaを使った方が楽
  - **chr, pos**, ref, alt
    - af_EUR
    - beta_EUR
    - se_EUR
    - **neglog10_pval_EUR**
    - low_confidence_EUR
  - **chr, pos**, ref, alt,
    - af_ (meta or each_of_6_pops)
    - beta_ <meta or each_of_6_pops>
    - se_ <meta or each_of_6_pops>
    - **neglog10_pval_** <meta or each_of_6_pops>
    - neglog10_pval_heterogeneity,
    - low_confidence_ <each_of_6_pops>
  - **chr, pos**, ref, alt,
    - af_ <meta_hq or meta or each_of_6_pops>
    - beta_ <meta_hq or meta or each_of_6_pops>
    - se_ <meta_hq or meta or each_of_6_pops>
    - **neglog10_pval_** <meta or each_of_6_pops>
    - neglog10_pval_ <heterogeneity_hq or heterogeneity
    - low_confidence_ <each_of_6_pops>

**finngen**
- 参考：https://www.finngen.fi/en/access_results
- ブラウザー：https://r11.finngen.fi/降しか無くて古いのはアクセスできない、左下のDownload Tableからcase数などの一覧が手に入る
- sumsta
  - finngen_xxx.gz
  - **#chrom**, **pos**, ref, alt, **rsids**, nearest_genes, **pval**, mlogp, beta, sebeta, af_alt, af_alt_cases, af_alt_controls

**Yang lab fastGWA**
- about: https://yanglab.westlake.edu.cn/data/ukb_fastgwa/imp_binary/about
- tutorial: https://yanglab.westlake.edu.cn/data/ukb_fastgwa/tutorial.html
- detail: https://yanglab.westlake.edu.cn/software/gcta/index.html#DataResource
  - sumstaの列の解説
  - Summary table: UKB_impute_v1.1.csv, UKB_binary_v1.11.csv, etc.
    - case数、control数  
- sumsta
  - <xxx>.v1.1.fastGWA.gz
    - **CHR**, **SNP**, **POS**, A1, A2, N, AF1, BETA, SE, **P**
  - <xxx>.v1.0.fastGWA.gz
    - **CHR**, **SNP**, **POS**, A1, A2, N, AF1, T, SE_T, P_noSPA, BETA, SE, **P**, CONVERGE
