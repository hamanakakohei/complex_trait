#!/usr/bin/env Rscript
#
# eQTL-eQTL coloc
#
# 出力：
# coloc.susie結果
#
# to do:
# - coloc.abfを引数で指定して使えるようにする
# - only_plot引数は廃止して、自動で毎回reginoal plotが出るようにする？
library(coloc)
library(arrow)
library(bigsnpr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--variant_id", help = "coloc疑いのセット、snp (phen_id1 or 2 のリード(indep)？) 例: chr22:17802533:T:C", required = TRUE)
parser$add_argument("--phen_id1", help = "coloc疑いのセット、gene1 例: Group224", required = TRUE)
parser$add_argument("--phen_id2", help = "coloc疑いのセット、gene2 例: ENSG00000015475.19", required = TRUE)
parser$add_argument("--assoc_prefix", help = "tensorqtl nominal結果の prefix（<prefix>chr??.parquet）（susie的にはphenファイルで良いのだが、coloc.susie的にはこっちのほうが早そう、、）", required = TRUE)
parser$add_argument("--geno_prefix", help = "bigsnpr用に変換した plink1 bed の prefix（例: path/to/genotypes.）、(tensorqtl時に合わせて、pgen/pvar/psam形式を受けたかったが、、）", required = TRUE)
parser$add_argument("--phen", help = "tensorqtl nominalで使ったphenotypeファイル、解析されたサンプルリストを作るためのみに使う", required = TRUE)
parser$add_argument("--cov", help = "tensorqtl nominalで使ったcovファイル、in-sample LDをregress outする", required = TRUE)
parser$add_argument("--window", help = "tensorqtl nominal sumsta（qtl1, qtl2）や続くやDをvariant_idを中心にしぼる", type = "integer", default = 1000000)
parser$add_argument("--coloc_type", help = "susie or abf", default = "susie")
parser$add_argument("--adjust_ld_by_cov", help = "covariteでadjustしたgenotypeでLDを計算する。", action = "store_true")
parser$add_argument("--out_prefix", help = ".is_in_cs95ファイルはcs列が-1ならcs95にない。とregional plot画像。とcoloc_susie結果のファイル。")
parser$add_argument("--only_plot", help = "ウィンドウを決めるためにプロットのみする。", action = "store_true")
args <- parser$parse_args()


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

make_D <- function(qtl, LD, N) {
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

make_D_wrap <- function(qtl, G, N, Z = NULL, adjust_covariates = TRUE){
  if (args$only_plot) {
    return( make_D(qtl, NULL, N) )
  }

  G_imp <- impute_mean(G)
  if (adjust_covariates) {
    G_use <- remove.covariate.effects(G_imp, Z, G_imp)$X
  }else{
    G_use <- G_imp
  }
  # 行名・列名の順番を保証
  LD <- cor(G_use)
  rownames(LD) <- colnames(G_use)
  colnames(LD) <- colnames(G_use)

  make_D(qtl, LD, N)
}

has_cs <- function(S) {
  !is.null(summary(S)$cs)
}

run_susie_with_fallback <- function(D1, D2, coverages = seq(0.95, 0.05, by = -0.05)) {
  for (cov in coverages) {
    message(sprintf("Trying coverage = %.2f", cov))
    S1 <- runsusie(D1, coverage = cov)
    S2 <- runsusie(D2, coverage = cov)

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


# 1. covの処理、転置してsamples × covariatesとするのとサンプルリストを作る
cov <- read.table(
  args$cov,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

cov_samples <- colnames(cov)
Z <- t(as.matrix(cov))


# 2. phenの処理、サンプルリストを作る
pheno <- read.table(
  args$phen,
  header = TRUE,
  comment.char = "",
  sep = "\t",
  check.names = FALSE
)

pheno_samples <- colnames(pheno)[5:ncol(pheno)]


# 3. tensorqtl nominal結果の処理、phen_id1/2のみのsumstaにする
chr <- sub(":.*", "", args$variant_id)
nominal_file <- paste0(args$assoc_prefix, chr, ".parquet")

lead_pos <- as.integer(sub(".*:(\\d+):.*", "\\1", args$variant_id))

nom <- read_parquet(nominal_file)
nom$pos <- as.integer(sub(".*:(\\d+):.*", "\\1", nom$variant_id))
nom <- nom[nom$pos >= lead_pos - args$window & nom$pos <= lead_pos + args$window, ]

qtl1 <- nom[nom$phenotype_id == args$phen_id1, ]
qtl2 <- nom[nom$phenotype_id == args$phen_id2, ]


# 4. genoの処理、plink bed → bigSNPと変換してgenotype matrixを得つつ、サンプルリストを作る
# bigsnprの中間ファイルをtmpdirに保存して最後に消す
tmpdir <- tempfile("bigsnpr_")
dir.create(tmpdir)

on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

backingfile <- file.path(tmpdir, "geno")

plink_bed <- paste0(args$geno_prefix, chr, ".bed")
snp_readBed(plink_bed, backingfile = backingfile)
obj <- snp_attach(paste0(backingfile, ".rds"))
G_all <- obj$genotypes
map   <- obj$map
fam   <- obj$fam
geno_samples <- fam$sample.ID


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


# 6. genoからcovをregress out --> in-sample LD計算 --> D作る
N = length(common_samples)
ind.col1 <- match(qtl1$variant_id, map$marker.ID)
G_mat1 <- G_mat[, ind.col1]
D1 <- make_D_wrap(qtl1, G_mat1, N, Z, adjust_covariates = COV_ADJUSTMENT)

ind.col2 <- match(qtl2$variant_id, map$marker.ID)
G_mat2 <- G_mat[, ind.col2]
D2 <- make_D_wrap(qtl2, G_mat2, N, Z, adjust_covariates = COV_ADJUSTMENT)


# 7. regional plot描いて、windowが適切か考える
if (args$only_plot) {
  xlim_common <- c(lead_pos - args$window / 3, lead_pos + args$window / 3)
  xticks <- seq(xlim_common[1], xlim_common[2], by = 2e5)

  png(paste0(args$out_prefix, ".png"), width = 3600, height = 1000, res = 120)
  par(mfrow = c(2, 1))
  plot_dataset(D1, main = args$variant_id, xlim = xlim_common, xaxt = "n")
  axis(1, at = xticks)
  plot_dataset(D2, main = "D2",            xlim = xlim_common, xaxt = "n")
  axis(1, at = xticks)
  dev.off()

  quit(status = 0)
}


# 8. susie
if(args$coloc_type == "susie"){

res_susie <- run_susie_with_fallback(D1, D2)
S1 <- res_susie$S1
S2 <- res_susie$S2
used_cov <- res_susie$coverage

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
saveRDS(used_cov, paste0(args$out_prefix, ".used_cov.rds"))


# 9. coloc
coloc_res <- coloc.susie(S1, S2)
saveRDS(coloc_res, paste0(args$out_prefix, ".coloc_res.rds"))

# パターン0：S2にCSないとき
if (!has_cs(S2)){

  # でもそもそもS1でVARIANT_IDがCSに含まれているかなぞ
  warning("No CS exists in D2 for the same coerage as D1")

} else if (used_cov == 0.95) {
  # パターン１：95%CSあるとき
  res <- coloc_res$summary

  ## オプション：coloc hitがターゲットSNPの場合のみにしたいとき
  ##res <- res[res$hit1 == args$variant_id, , drop = FALSE]
  #res <- res[res$hit1 == VARIANT_ID, , drop = FALSE]
  #
  #if (nrow(res) == 0) {
  #  warning("Target variant not  found in coloc.susie summary")
  #}

  res$gene1 <- args$phen_id1
  res$gene2 <- args$phen_id2
  res$cs    <- used_cov
  write.table(
    res,
    file = paste0(args$out_prefix, ".coloc_susie_res.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
} else if (used_cov != 0.95) {
  # パターン２：95%CSないとき
  res <- coloc_res$results

  ## オプション：coloc hitがターゲットSNPの場合のみにしたいとき
  ##res <- res[res$snp == args$variant_id, , drop = FALSE]
  #res <- res[res$snp == VARIANT_ID, , drop = FALSE]
  #
  #if (nrow(res) == 0) {
  #  stop("Target variant not found in coloc.susie results")
  #}

  # summary にしか無い列を NA で埋める
  res$nsnps            <- NA
  res$hit1             <- args$variant_id
  res$hit2             <- args$variant_id
  res$PP.H0.abf        <- NA
  res$PP.H1.abf        <- NA
  res$PP.H2.abf        <- NA
  res$PP.H3.abf        <- NA
  res$idx1             <- NA
  res$idx2             <- NA
  res$gene1            <- args$phen_id1
  res$gene2            <- args$phen_id2
  res$cs               <- used_cov
  write.table(
    res,
    file = paste0(args$out_prefix, ".coloc_susie_res.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
