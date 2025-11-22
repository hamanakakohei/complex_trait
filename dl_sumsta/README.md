gwas catalog
ieu openGWAS

neale lab ukbb
- xxx_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
  - variant, minor_allele, minor_AF, low_confidence_variant, n_complete_samples, AC, ytx, beta, se, tstat, pval
- biomarkers-xxx-both_sexes-irnt.tsv.bgz 
  - chr, pos, ref, alt
    - af_<meta_hq or meta or each of 6 pops>
    - beta_<meta_hq or meta or each of 6 pops>
    - se_<meta_hq or meta or each of 6 pops>
    - neglog10_pval_<meta_hq or meta or each of 6 pops>
    - neglog10_pval_heterogeneity_hq, neglog10_pval_heterogeneity
    - low_confidence_<each_of_6_pops>
- categorical-xxx-both_sexes-xxx.tsv.bgz
  - chr, pos, ref, alt
    - af_cases_<meta or each_of_6_pops>
    - af_controls_<meta or each_of_6_pops>
    - beta_<meta or each_of_6_pops>
    - se_<meta or each_of_6_pops>
    - neglog10_pval_<meta or each_of_6_pops>
    - neglog10_pval_heterogeneity
    - low_confidence_<each_of_6_pops>
- continuous-xxx--both_sexes-<"" or "auto_medadj" or "combined_medadj">_<irnt or raw>.tsv.bgz
  - chr,pos,ref,alt,af_EUR,beta_EUR,se_EUR,neglog10_pval_EUR,low_confidence_EUR
  - chr,pos,ref,alt,
    - af_<meta or each_of_6_pops>
    - beta_<meta or each_of_6_pops>,
    - se_<meta or each_of_6_pops>,
    - neglog10_pval_<meta or each_of_6_pops>,
    - neglog10_pval_heterogeneity,
    - low_confidence_<each_of_6_pops>
  - chr,pos,ref,alt,
    - af_<meta_hq or meta or each of 6 pops>
    - beta_<meta_hq or meta or each of 6 pops>
    - se_<meta_hq or meta or each of 6 pops>
    - neglog10_pval_<meta or each_of_6_pops>
    - neglog10_pval_heterogeneity_hq, neglog10_pval_heterogeneity
    - low_confidence_<each_of_6_pops>

finngen_xxx.gz
- #chrom, pos, ref, alt, rsids, nearest_genes, pval, mlogp, beta, sebeta, af_alt, af_alt_cases, af_alt_controls

Yang lab
- <xxx>.v1.1.fastGWA.gz
  - CHR, SNP, POS, A1, A2, N, AF1, BETA, SE, P
- <xxx>.v1.0.fastGWA.gz
  - CHR, SNP, POS, A1, A2, N, AF1, T, SE_T, P_noSPA, BETA, SE, P, CONVERGE
