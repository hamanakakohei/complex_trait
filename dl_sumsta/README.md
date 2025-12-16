to do
- gwas catalog
  - ダウンロードをapiでする
  - urlが間違っていてもgcst idから自動で出来るようにする
- 他のdbも既製のスクリプトがないのか？？

**GeneATLAS**
- trait key表：http://geneatlas.roslin.ed.ac.uk/traits-table/
- ダウンロード用スクリプトはhttp://geneatlas.roslin.ed.ac.uk/downloads/で適当なtraitを入れると汎用的なものが得られてそれを使い回す
- 列（参照：http://www.dissect.ed.ac.uk/documentation-gwas/
-   SNP ALLELE iscores NBETA-xxx NSE-xxx PV-xxx

ieu openGWAS

**gwas catalog**
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
- 参照：https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=227859291#gid=227859291
- xxx_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
  - **variant**, minor_allele, minor_AF, low_confidence_variant, n_complete_samples, AC, ytx, beta, se, tstat, **pval**

**pan-ukbb**
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
- finngen_xxx.gz
  - **#chrom**, **pos**, ref, alt, **rsids**, nearest_genes, **pval**, mlogp, beta, sebeta, af_alt, af_alt_cases, af_alt_controls

**Yang lab**
- <xxx>.v1.1.fastGWA.gz
  - **CHR**, **SNP**, **POS**, A1, A2, N, AF1, BETA, SE, **P**
- <xxx>.v1.0.fastGWA.gz
  - **CHR**, **SNP**, **POS**, A1, A2, N, AF1, T, SE_T, P_noSPA, BETA, SE, **P**, CONVERGE
