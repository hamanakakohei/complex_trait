#!/usr/bin/env Rscript

library(tidyverse)
library(argparser)
library(exact2x2)


arg_parser("分割表からフォレプロを描く、すでにORなど計算されている版も作りたい") %>%
  add_argument( "--contingency_table", help = "bin列、data列ごとのn数") %>%
  add_argument( "--effect_col", help = "ORの分子の列") %>%
  add_argument( "--noneffect_col", help = "ORの分母の列") %>%
  add_argument( "--ref_bin", help = "ORの基準bin") %>%
  add_argument( "--bin_order", nargs="+", help = "上のtable内ののbin列のプロット順、最初にref_binをもってくること") %>%
  add_argument( "--data_order", nargs="+", help = "上のtable内ののdata列のプロット順") %>%
  add_argument( "--out", help="") %>%
  parse_args() -> argv

argv$bin_order <- strsplit(argv$bin_order, ",")[[1]]
argv$data_order <- strsplit(argv$data_order, ",")[[1]]

dt = read_tsv(argv$contingency_table)


# 1：ORなどを計算する
plot_bins = setdiff(argv$bin_order, argv$ref_bin)

all_res = lapply( unique(dt$data), function( DATA ){
  lapply( plot_bins, function( PLOT_BIN ){

      table.2x2 = dt %>% 
        mutate( bin = factor( bin, levels = argv$bin_order ) ) %>%
        arrange( bin ) %>%
        filter( data == DATA ) %>% 
        filter( bin %in% c( argv$ref_bin, PLOT_BIN ) ) %>% 
        select( argv$noneffect_col, argv$effect_col ) %>%
        as.matrix
      print(table.2x2)
      if (nrow(table.2x2) != 2 || ncol(table.2x2) != 2)
        stop("Error: 分割表が変です")
      
      res = fisher.exact( table.2x2 )
    
      tibble(
        data = DATA,
        bin = PLOT_BIN,
        or = as.numeric(res$estimate),
        dn = as.numeric(res$conf.int[1]),
        up = as.numeric(res$conf.int[2]),
        p = res$p.value
      )
    }) %>%
    do.call( rbind, . )

  }) %>%
  do.call( rbind, . )


# 2：プロット
p = ggplot( all_res, aes( x = data, y = or, ymin = dn, ymax = up, group = bin, color = bin ) ) +
  geom_pointrange( position = position_dodge( width = 0.6 ) ) +  # aes( col = bin ), lwd = 0.8
  geom_errorbar( width = 0.3, position = position_dodge( width = 0.6 )) + # aes( ymin = dn, ymax = up, col = bin ), cex = 1
  geom_hline( yintercept = 1, linetype = "dashed" ) +
  coord_flip() +
  scale_x_discrete( limits = argv$data_order ) +
  theme_minimal()
  #theme(
  #  strip.text.y = element_text( hjust=0, vjust = 1, angle=180, face="bold")
  #) +
  #scale_color_manual( values = c(
  #  "1" = rgb(0,0,0,0.1),
  #  "2" = rgb(0,0,0,0.25),
  #  "3" = rgb(0,0,0,1.0)
  #) ) +

ggsave( argv$out, p ) #, height=60, units="cm", dpi=600 )

