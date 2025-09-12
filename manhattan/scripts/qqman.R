#!/usr/bin/env Rscript
library(qqman)

# ダミーデータ作成
set.seed(123)
n_genes <- 50000
genes <- paste0("Gene", 1:n_genes)
chroms <- sample(1:22, n_genes, replace=TRUE)
pos <- sample(1:1e7, n_genes, replace=TRUE)

# P値: 多くは適当に大きめ、少数だけ強いシグナルを混ぜる
# 基本的には有意でない
pvals <- runif(n_genes, min=1e-5, max=1) 
# 散発的なシグナル
n_genes_sig = 50
signal_genes <- sample(1:n_genes, n_genes_sig)    
# そのままで均等版
#pvals[signal_genes] <- runif(n_genes_sig, min=1e-10, max=1e-5)
# -log10して均等版
log10p <- runif(n_genes_sig, min=5, max=10) 
pvals[signal_genes] <- 10^(-log10p)      

# データフレームをqqman形式に
df <- data.frame(
  SNP = genes,
  CHR = chroms,
  BP  = pos,
  P   = pvals
)

# Manhattan plot
pdf("qqman.pdf", width=4, height=3)
manhattan(df,
          main="Gene-based Manhattan plot (dummy)",
          ylim=c(0, 10),
          genomewideline=FALSE, # -log10(5e-8), 
          suggestiveline=FALSE, # -log10(1e-5), 
          cex=0.6, col=c("steelblue","darkorange"))
dev.off()
