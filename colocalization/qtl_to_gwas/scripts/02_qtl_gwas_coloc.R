#!/usr/bin/env Rscript
#
# eQTL-GWAS coloc
#
# 出力：
# coloc.susie結果
library(coloc)
library(arrow)
library(bigsnpr)
library(argparse)
library(ieugwasr)

parser <- ArgumentParser()
parser$add_argument("--variant_id", help = "coloc疑いのセット、snp (phen_id1 or 2 のリード(indep)？) 例: chr22:17802533:T:C", required = TRUE)
parser$add_argument("--egene", help = "coloc疑いのeGene。例: ENSG00000015475.19", required = TRUE)
parser$add_argument("--eqtl_assoc_prefix", help = "tensorqtl nominal結果の prefix（<prefix>chr??.parquet）（susie的にはphenファイルで良いのだが、coloc.susie的にはこっちのほうが早そう、、）", required = TRUE)
parser$add_argument("--eqtl_geno_prefix", help = "bigsnpr用に変換した plink1 bed の prefix（例: path/to/genotypes.）、(tensorqtl時に合わせて、pgen/pvar/psam形式を受けたかったが、、）", required = TRUE)
parser$add_argument("--eqtl_phen", help = "tensorqtl nominalで使ったphenotypeファイル、解析されたサンプルリストを作るためのみに使う", required = TRUE)
parser$add_argument("--eqtl_cov", help = "tensorqtl nominalで使ったcovファイル、in-sample LDをregress outする", required = TRUE)
parser$add_argument("--window", help = "tensorqtl nominal sumsta（qtl1, qtl2）や続くやDをvariant_idを中心にしぼる", type = "integer", default = 1000000)
parser$add_argument("--gwas_sumsta", help = "", required = TRUE)
parser$add_argument("--gwas_type", help = "cc or quant", required = TRUE)
parser$add_argument("--gwas_ld_pop", help = "", required = TRUE)
parser$add_argument("--gwas_n", type = "integer", required = TRUE)
parser$add_argument("--gwas_n_case", type = "integer", help = "gwas_typeがccなら必要")
parser$add_argument("--plink_bin", default = "/home/khamanaka/.local/bin/plink_downloaded/plink1.9")
parser$add_argument("--ld_ref_bfile_prefix", default = "/home/khamanaka/resource/1kg/plink_homepage/all_hg38_autosome_")
#parser$add_argument("--ld_ref_bfile_prefix", default = "/home/khamanaka/resource/1kg/1kg.v3__hg19/")
parser$add_argument("--out_prefix", help = ".is_in_cs95ファイルはcs列が-1ならcs95にない。とregional plot画像。とcoloc_susie結果のファイル。")
parser$add_argument("--only_plot", help = "ウィンドウを決めるためにプロットのみする。", action = "store_true")
args <- parser$parse_args()

VARIANT_ID = args$variant_id
EGENE = args$egene
EQTL_ASSOC_PREFIX = args$eqtl_assoc_prefix
EQTL_GENO_PREFIX = args$eqtl_geno_prefix
EQTL_PHEN = args$eqtl_phen
EQTL_COV = args$eqtl_cov
WINDOW = args$window
GWAS_SUMSTA = args$gwas_sumsta
GWAS_TYPE = args$gwas_type
GWAS_LD_POP = args$gwas_ld_pop
GWAS_N = args$gwas_n
GWAS_N_CASE = args$gwas_n_case
ONLY_PLOT = args$only_plot
PLINK_BIN = args$plink_bin
LD_REF_BFILE_PREFIX = args$ld_ref_bfile_prefix

print(VARIANT_ID)
print(EGENE)
print(EQTL_ASSOC_PREFIX)
print(EQTL_GENO_PREFIX)
print(EQTL_PHEN)
print(EQTL_COV)
print(WINDOW)
print(GWAS_SUMSTA)
print(GWAS_TYPE)
print(GWAS_LD_POP)
print(GWAS_N)
print(GWAS_N_CASE)
print(ONLY_PLOT)
print(PLINK_BIN)
print(LD_REF_BFILE_PREFIX)

## 5
#VARIANT_ID = 'chr11:111378732:A:T'
#EGENE = 'Group165'
#
#EQTL_ASSOC_PREFIX = '~/gm_rnaseq/gtex-pipeline_do2/qtl_e/cromwell-executions/tensorqtl_cis_nominal_workflow/03193177-34d0-4925-b899-314cb653c348/call-tensorqtl_cis_nominal/execution/rnaseq430.nominal.cis_qtl_pairs.'
#EQTL_GENO_PREFIX = '~/gm/plink_bed/plink1_format/results/ALL.correctRefAlt.norm.'
#EQTL_PHEN = '~/gm_rnaseq/gtex-pipeline_do2/qtl_e/results/04/rnaseq430.expression.bed.gz'
#EQTL_COV = '~/gm_rnaseq/gtex-pipeline_do2/qtl_e/results/04/rnaseq430.combined_covariates.txt'
#WINDOW = 300000
#GWAS_SUMSTA = met-d-XL_HDL_C_pct.vcf.gz
#GWAS_LD_POP =
#GWAS_N =
#GWAS_N_CASE =
#ONLY_PLOT = FALSE


# 0. 関数
impute_mean <- function(G) {
  for (j in seq_len(ncol(G))) {
    idx <- is.na(G[, j])
    if (any(idx)) {
      G[idx, j] <- mean(G[, j], na.rm = TRUE)
    }
  }
  G
}

remove.covariate.effects <- function (X, Z, y) {
  if (any(Z[,1] != 1)) Z <- cbind(1, Z)
  A   <- crossprod(Z)
  SZy <- solve(A, crossprod(Z, y))
  SZX <- solve(A, crossprod(Z, X))
  y   <- y - Z %*% SZy
  X   <- X - Z %*% SZX
  list(X = X, y = y)
}

make_D_4_qtl <- function(qtl, LD, N) {
  pos <- as.integer(sub(".*:(\\d+):.*", "\\1", qtl$variant_id))

  list(
    beta    = qtl$slope,
    varbeta = qtl$slope_se^2,
    N       = N,
    type    = "quant",
    MAF     = qtl$af,
    LD      = LD,
    snp     = qtl$variant_id,
    position = pos
  )
}

make_D_4_qtl_wrap <- function(qtl, G, N){
  #if (args$only_plot) {
  if (ONLY_PLOT) {
    make_D_4_qtl(qtl, NULL, N)
  }else{
    G_imp <- impute_mean(G)
    G_resid <- remove.covariate.effects(G_imp, Z, G_imp)$X

    # 行名・列名の順番を保証
    LD <- cor(G_resid)
    rownames(LD) <- colnames(G_resid)
    colnames(LD) <- colnames(G_resid)

    make_D_4_qtl(qtl, LD, N)
  }
}

make_D_4_gwas <- function(gwas_sumsta, gwas_type, LD, N, N_CASE=NULL) {
  x <- list(
    beta     = gwas_sumsta$BETA,
    varbeta  = gwas_sumsta$VAR_BETA,
    N        = N,
    LD       = LD,
    snp      = gwas_sumsta$SNPID,
    position = gwas_sumsta$POS
  )

  if (gwas_type == "quant") {
    x$type <- "quant"
    if ("EAF" %in% colnames(gwas_sumsta)) {
      x$MAF <- gwas_sumsta$EAF
    }else{
      x$sdY <- 1
    }
  } else if (gwas_type == "cc") {
    x$type <- "cc"
    x$s    <- N_CASE / N
  }

  return(x)
}

make_D_4_gwas_wrap <- function(gwas_sumsta, gwas_type, LD, N, N_CASE=NULL){
  make_D_4_gwas(gwas_sumsta, gwas_type, LD, N, N_CASE)
}

has_cs <- function(S) {
  !is.null(summary(S)$cs)
}

run_susie_with_coverage_fallback <- function(D1, D2, coverages = seq(0.95, 0.05, by = -0.05)) {
  for (cov in coverages) {
    message(sprintf("Trying coverage = %.2f", cov))
    S1 <- runsusie(D1, coverage = cov, maxit=10000, repeat_until_convergence=FALSE)
    S2 <- runsusie(D2, coverage = cov, maxit=10000, repeat_until_convergence=FALSE)

    if (has_cs(S1) && has_cs(S2)) {
      message(sprintf("Success at coverage = %.2f", cov))
      return(list(
        S1 = S1,
        S2 = S2,
        coverage = cov
      ))
    }
  }

  # すべて失敗
  stop("No credible sets found for any coverage")
}

run_susie_with_abf_fallback <- function(D1, D2, cov=0.95) {
  susie_ok <- TRUE
  S1 <- S2 <- NULL

  tryCatch({
    S1 <- runsusie(D1, coverage = cov, maxit=10000, repeat_until_convergence=FALSE)
    S2 <- runsusie(D2, coverage = cov, maxit=10000, repeat_until_convergence=FALSE)
    for(i in 1:10){ trash <- summary(S1) }
    for(i in 1:10){ trash <- summary(S2) }
  }, error = function(e) {
    susie_ok <<- FALSE
  })

  if (susie_ok && has_cs(S1) && has_cs(S2)) {
    message("susie success")
    print("D1 susie res summary")
    print(summary(S1))
    print("D1 susie cs")
    print(S1$sets)
    print("D2 susie res summary")
    print(summary(S2))
    print("D2 susie cs")
    print(S2$sets)

    saveRDS(S1, paste0(args$out_prefix, ".S1.rds"))
    saveRDS(S2, paste0(args$out_prefix, ".S2.rds"))

    coloc_res <- coloc.susie(S1, S2)
    saveRDS(coloc_res, paste0(args$out_prefix, ".coloc_res.rds"))

    res <- coloc_res$summary
    #res$egene <- args$egene
    #res$gwas  <- args$gwas_sumsta
    res$egene <- EGENE
    res$gwas  <- GWAS_SUMSTA
    res$coloc_type  <- "susie"
    write.table(
      res,
      file = paste0(args$out_prefix, ".coloc_res.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  } else {
    message("abf selected")
    coloc_res <- coloc.abf(D1, D2)
    saveRDS(coloc_res, paste0(args$out_prefix, ".coloc_res.rds"))

    print(head((coloc_res$results[order(coloc_res$results$SNP.PP.H4, decreasing=T),])))

    res <- coloc_res$summary
    #res$egene <- args$egene
    #res$gwas  <- args$gwas_sumsta
    res$hit1  <- NA
    res$hit2  <- NA
    res$idx1  <- NA
    res$idx2  <- NA
    res$egene <- EGENE
    res$gwas  <- GWAS_SUMSTA
    res$coloc_type  <- "abf"
    write.table(
      res,
      file = paste0(args$out_prefix, ".coloc_res.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}


# 1. covの処理、転置してsamples × covariatesとするのとサンプルリストを作る
cov <- read.table(
  #args$eqtl_cov,
  EQTL_COV,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

cov_samples <- colnames(cov)
Z <- t(as.matrix(cov))
print("1111111111")


# 2. phenの処理、サンプルリストを作る
pheno <- read.table(
  #args$eqtl_phen,
  EQTL_PHEN,
  header = TRUE,
  comment.char = "",
  sep = "\t",
  check.names = FALSE
)

pheno_samples <- colnames(pheno)[5:ncol(pheno)]
print("222222222")


# 3. tensorqtl nominal結果の処理、phen_id1/2のみのsumstaにする
#chr <- sub(":.*", "", args$variant_id)
chr <- sub(":.*", "", VARIANT_ID)
#nominal_file <- paste0(args$eqtl_assoc_prefix, chr, ".parquet")
nominal_file <- paste0(EQTL_ASSOC_PREFIX, chr, ".parquet")

#lead_pos <- as.integer(sub(".*:(\\d+):.*", "\\1", args$variant_id))
lead_pos <- as.integer(sub(".*:(\\d+):.*", "\\1", VARIANT_ID))

nom <- read_parquet(nominal_file)
nom$pos <- as.integer(sub(".*:(\\d+):.*", "\\1", nom$variant_id))
#nom <- nom[nom$pos >= lead_pos - args$window & nom$pos <= lead_pos + args$window, ]
nom <- nom[nom$pos >= lead_pos - WINDOW & nom$pos <= lead_pos + WINDOW, ]

#qtl <- nom[nom$phenotype_id == args$egene, ]
qtl <- nom[nom$phenotype_id == EGENE, ]
print("333333333")


# 4. genoの処理、plink bed → bigSNPと変換してgenotype matrixを得つつ、サンプルリストを作る
# bigsnprの中間ファイルをtmpdirに保存して最後に消す
tmpdir <- tempfile("bigsnpr_")
dir.create(tmpdir)

on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

backingfile <- file.path(tmpdir, "geno")

#plink_bed <- paste0(args$eqtl_geno_prefix, chr, ".bed")
plink_bed <- paste0(EQTL_GENO_PREFIX, chr, ".bed")
snp_readBed(plink_bed, backingfile = backingfile)
obj <- snp_attach(paste0(backingfile, ".rds"))
G_all <- obj$genotypes
map   <- obj$map
fam   <- obj$fam
geno_samples <- fam$sample.ID
print("4444444444")


# 5. cov, phen, genoのサンプルリストからcommon samplesを得て、covとgenoをそれにそろえつつ、genoはmatrix化して列名行名を付ける
common_samples <- Reduce(
  intersect,
  list(cov_samples, pheno_samples, geno_samples)
)

ind.row <- match(common_samples, fam$sample.ID)

G <- G_all[ind.row, ]
G_mat <- as.matrix(G)
rownames(G_mat) <- common_samples
colnames(G_mat) <- map$marker.ID

Z <- Z[common_samples, , drop = FALSE]
print("55555555")


# 6. genoからcovをregress out --> in-sample LD計算 --> D作る手前まで
QTL_N = length(common_samples)
ind.col1 <- match(qtl$variant_id, map$marker.ID)
G_mat1 <- G_mat[, ind.col1]
print("66666666")


# 6-2. gwas側のD作る
# gwas sumstaをqtl variantに限る、じゃないとld計算つらい
gwas <- read.table(GWAS_SUMSTA, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
gwas <- gwas[!duplicated(gwas$SNPID), ] # GWASlabでやるべき、一部のデータでSNPが被っている、、他の列の値は微妙に違う、、
print(head(gwas))
print(sapply(gwas, class))
if ("EAF" %in% colnames(gwas)) {
  gwas <- gwas[gwas$EAF > 0 & gwas$EAF < 1, ] # GWASlabでやるべき
}
print(head(gwas))
gwas$SNPID <- paste0("chr", gwas$SNPID)
head(dim(gwas))
head(gwas)
cat("gwas n: ", nrow(gwas), "\n")
gwas <- gwas[gwas$SNPID %in% qtl$variant_id, ]
cat("qtl$variant_id: ", length(qtl$variant_id), "\n")
cat("gwas n (only qtl variant): ", nrow(gwas), "\n")
cat("gwas n (only qtl variant) examples: ", head(gwas$SNPID), "\n")
head(dim(gwas))
head(gwas)
print("777777777777")

# ld計算してld matの行 or 列の要素がすべて NA なら除き
ld_mat <- ld_matrix(
  gwas$rsID,
  plink_bin = PLINK_BIN,
  bfile = paste0(LD_REF_BFILE_PREFIX, GWAS_LD_POP)
)
cat("ld_mat n (only gwas variant): ", nrow(ld_mat), "\n")
print(dim(ld_mat))
print(ld_mat[1:5,1:5])
all_na <- apply(ld_mat, 1, function(x) all(is.na(x)))
ld_mat <- ld_mat[!all_na, !all_na]
cat("ld_mat n (after removing NA): ", nrow(ld_mat), "\n")
print(dim(ld_mat))
print(ld_mat[1:5,1:5])
if (any(is.na(ld_mat))) {
  stop("LD matrix still contains NA after removing all-NA SNPs")
}
print("8888888")

# ld計算のrsID, ref, altがgwas側の想定と合っていることを確認し
# （合っていないならld rをフリップしたいが未実装）
# ld matのsnp順をgwas側に合わせる
ld_keys <- do.call(rbind, strsplit(rownames(ld_mat), "_"))
ld_snp_df <- data.frame(
  rsID = ld_keys[,1],
  REF  = ld_keys[,3], # ここ注意! 2でなくて3
  ALT  = ld_keys[,2],
  stringsAsFactors = FALSE
)
print(dim(ld_snp_df))
print(head(ld_snp_df))

merged <- merge(gwas, ld_snp_df, by.x = c("rsID", "NEA", "EA"), by.y = c("rsID", "REF", "ALT"), all = FALSE)
cat("merged n (not flipped allele btw gwas & ld_mat)", nrow(merged), "\n")
merged$ld_row <- match(
  paste(merged$rsID, merged$EA, merged$NEA, sep = "_"),
  paste(ld_snp_df$rsID, ld_snp_df$ALT, ld_snp_df$REF, sep = "_")
)

if (any(is.na(merged$ld_row))) {
  stop("Some merged SNPs not found in LD matrix")
}
print("9999999999")

gwas <- merged
print(dim(gwas))
print(head(gwas))
ld_mat <- ld_mat[merged$ld_row, merged$ld_row]
print(dim(ld_mat))
print(ld_mat[1:5,1:5])
colnames(ld_mat) <- gwas$SNPID
rownames(ld_mat) <- gwas$SNPID
print(dim(ld_mat))
print(ld_mat[1:5,1:5])
print(gwas$SNPID[duplicated(gwas$SNPID)])
print("00000000000")


# 6-3. リードが含まれるか？D1とD2のsnpの被り数は？
# shared snpのみにしないといけないかは確認していないが、snpの順番をD1とD2ではそろえないといけないようだ
shared_snp = intersect(qtl$variant_id, gwas$SNPID)

D1_raw <- make_D_4_qtl_wrap(qtl, G_mat1, QTL_N)
D2_raw <- make_D_4_gwas_wrap(gwas, GWAS_TYPE, ld_mat, GWAS_N, GWAS_N_CASE)

saveRDS(D1_raw, paste0(args$out_prefix, ".D1_raw.rds"))
saveRDS(D2_raw, paste0(args$out_prefix, ".D2_raw.rds"))

D1 <- make_D_4_qtl_wrap(
        qtl[match(shared_snp, qtl$variant_id), ],
        G_mat1[, match(shared_snp, colnames(G_mat1))],
        QTL_N)
D2 <- make_D_4_gwas_wrap(
        gwas[match(shared_snp, gwas$SNPID), ],
        GWAS_TYPE,
        ld_mat[match(shared_snp, rownames(ld_mat)),
        match(shared_snp, colnames(ld_mat))],
        GWAS_N,
        GWAS_N_CASE)

cat(
  "D1_raw$snp:", length(D1_raw$snp), "\n",
  "D2_raw$snp:", length(D2_raw$snp), "\n",
  "Common:",     length(shared_snp), "\n",
  "D1$snp:",     length(D1$snp), "\n",
  "D2$snp:",     length(D2$snp), "\n"
)

if (!(VARIANT_ID %in% qtl$variant_id && VARIANT_ID %in% gwas$SNPID)) {
  stop("VARIANT_ID is missing in D1$snp or D2$snp")
}

saveRDS(D1, paste0(args$out_prefix, ".D1.rds"))
saveRDS(D2, paste0(args$out_prefix, ".D2.rds"))


# 7. regional plot描いて、windowが適切か考える
#if (args$only_plot) {
if (ONLY_PLOT) {
  xlim_common <- c(lead_pos - args$window / 3, lead_pos + args$window / 3)
  xticks <- seq(xlim_common[1], xlim_common[2], by = 2e5)

  png(paste0(args$out_prefix, ".png"), width = 3600, height = 1000, res = 120)
  par(mfrow = c(4, 1))
  plot_dataset(D1_raw, main = args$variant_id, xlim = xlim_common, xaxt = "n")
  axis(1, at = xticks)
  plot_dataset(D2_raw, main = "D2_raw",        xlim = xlim_common, xaxt = "n")
  axis(1, at = xticks)
  plot_dataset(D1,     main = "D1_shared",     xlim = xlim_common, xaxt = "n")
  axis(1, at = xticks)
  plot_dataset(D2,     main = "D2_shared",     xlim = xlim_common, xaxt = "n")
  axis(1, at = xticks)
  dev.off()

  quit(status = 0)
}


# 8. coloc
run_susie_with_abf_fallback(D1, D2)



#gwas2 <- gwas[gwas$rsID %in% ld_info$rsID, ]
#print(head(gwas))
#print(head(ld_info))
#print(head(gwas[gwas$rsID=="rs4911130",]))
#print(head(ld_info[ld_info$rsID=="rs4911130",]))
#print(head(gwas2))
#print(head(gwas2[gwas2$rsID=="rs4911130",]))
#
## gwas sumstaとLD matでアレルが逆になっていたらLD mat rの正負を逆にする
#flip <- logical(nrow(ld_info2))
#
#for (i in seq_len(nrow(ld_info2))) {
#
#  g_ea  <- gwas2$EA[i]
#  g_nea <- gwas2$NEA[i]
#
#  ld_ref <- ld_info$REF[i]
#  ld_alt <- ld_info$ALT[i]
#
#  if (g_ea == ld_alt && g_nea == ld_ref) {
#    flip[i] <- FALSE
#  } else if (g_ea == ld_ref && g_nea == ld_alt) {
#    flip[i] <- TRUE
#  } else {
#    stop(
#      sprintf(
#        "Allele mismatch at %s (GWAS %s/%s vs LD %s/%s)",
#        ld_info$rsID[i],
#        g_nea, g_ea,
#        ld_ref, ld_alt
#      )
#    )
#  }
#}
#
#flip_idx <- which(flip)
#
#if (length(flip_idx) > 0) {
#  ld_mat[flip_idx, ] <- -ld_mat[flip_idx, ]
#  ld_mat[, flip_idx] <- -ld_mat[, flip_idx]
#}
