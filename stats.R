#!/usr/bin/env Rscript

library(vroom)
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
library(tidyr)


if (!exists("snakemake")) {
  setClass("snakemake_fake", representation(
    params = "list", input = "list", output = "list",
    log = "list", wildcards = "list"))
  snakemake <- new("snakemake_fake",
    input = list("output/aj-enriched_v3.1.1-b38_chrall.tsv"),
    params = list(),
    log = list(),
    output = list("output/aj-enriched_v3.1.1-b38_an-ge20.stats.tsv"),
    wildcards = list(
      an = "20",
      build = "v3.1.1-b38"
    )
  )
}


extr_build <- function(x) strsplit(x, "-")[[1]]
build <- snakemake@wildcards[["build"]] %>% extr_build
build_gnomad <- build[1]
build_genome <- build[2]
an_aj <- as.integer(snakemake@wildcards[["an"]])

cs <-  cols(
  .default = col_character(),
  CHROM = col_number(),
  POS = col_integer(),
  AF_nfe = col_double(),
  AF_asj = col_double(),
  AN_asj = col_double()
)

stats <- vroom(snakemake@input[[1]], col_types = cs) %>%
  rename(SNP = ID, CHR = CHROM, BP = POS) %>%
  filter(AN_asj >= an_aj) %>%
  mutate(MAF_nfe = ifelse(AF_nfe > 0.5, 1 - AF_nfe, AF_nfe),
         MAF_asj = ifelse(AF_asj > 0.5, 1 - AF_asj, AF_asj),
         CHR = as.integer(CHR))

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

do_stats <- function(df) {
  summarise(df, N = n(), mean_nfe = mean(MAF_nfe),
            mean_asj = mean(MAF_asj),
            mode_nfe = getmode(MAF_nfe))
}

stats %>%
  group_by(CHR) %>%
  do_stats %>%
  mutate(CHR = as.character(CHR)) %>%
  bind_rows(mutate(do_stats(stats), CHR = "all")) %>%
  mutate(mean_nfe = format(mean_nfe, scientific = T, digits = 3),
         mode_nfe = ifelse(mode_nfe == 0, 0,
                          format(mode_nfe, scientific = T, digits = 3)),
         mean_asj = format(mean_asj, scientific = F, digits = 3)) %>%
  vroom_write(snakemake@output[[1]])
