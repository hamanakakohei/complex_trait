#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(locuscomparer)
  library(ggplot2)
})

source("scripts/locuscompareR_patch.R")
assignInNamespace("get_position", get_position, ns = "locuscomparer")
assignInNamespace("retrieve_LD", retrieve_LD, ns = "locuscomparer")


parser <- ArgumentParser(description = "")
parser$add_argument("--sumsta1",  required=TRUE, help="rsid, pvalという2列を想定")
parser$add_argument("--sumsta2",  required=TRUE, help="rsid, pvalという2列を想定")
parser$add_argument("--title1", default="1st_data")
parser$add_argument("--title2", default="2nd_data")
parser$add_argument("--rsid",   required=TRUE, help="lead")
parser$add_argument("--pop",   required=TRUE, choices=c("EUR", "EAS", "AFR", "SAS", "AMR", "GBR"))
parser$add_argument("--assembly",   required=TRUE, choices=c("hg19", "hg38"))
parser$add_argument("--out",   required=TRUE)
args <- parser$parse_args()


p <- locuscompare(
  in_fn1 = args$sumsta1,
  in_fn2 = args$sumsta2,
  title  = args$title1,
  title2 = args$title2,
  snp = args$rsid,
  population = args$pop,
  genome = args$assembly
)

ggsave(filename = args$out, plot = p, width = 8, height = 4, dpi = 300)
