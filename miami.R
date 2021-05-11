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
    output = list("output/aj-enriched_v3.1.1-b38_an-ge20.png"),
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

ss <- stats %>%
  pivot_longer(starts_with("MAF_"), names_to = "ancestry",
               names_prefix = "MAF_", values_to = "MAF") %>%
  mutate(MAF_miami = ifelse(ancestry == "asj", MAF, -MAF)) %>%
  select(CHR, POS = BP, MAF_miami)

don <- ss %>%
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len = max(POS), .groups = "drop") %>%
  mutate(chr_len = as.numeric(chr_len)) %>%
  # Calculate cumulative POS of each CHR
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(ss, ., by = c("CHR" = "CHR")) %>%
  # Add a cumulative POS of each SNP
  arrange(CHR, POS) %>%
  #mutate(BPcum = 1:nrow(dat)) ## No gaps in chromosome
  mutate(BPcum = POS + tot) ## Gaps in chromsome

message("Generating Miami plot")

plot_miami <- function(sumstats) {
  axisdf <- sumstats %>%
    group_by(CHR) %>%
    dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2, .groups = "drop")

  sumstats %>%
    mutate(chr.col = ifelse(as.numeric(as.character(CHR)) %% 2 == 1,
                            "odd", "even")) %>%
    ggplot() +
      geom_point(aes(x = BPcum, y = MAF_miami, colour = chr.col), size = 0.25) +
      coord_cartesian(ylim = c(-0.11, 0.51)) +
      scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
      scale_y_continuous(labels = abs, breaks = seq(-0.1, 0.5, 0.1)) +
      theme_bw() +
      theme(
        text = element_text(size = 10),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        ) +
      labs(y = "(non-Finnish European)    <- MAF ->    (Ashkenazi Jewish)",
           title = "Miami plot of MAF",
           subtitle = sprintf("gnomad %s (%s); Minumum AJ Allele Number = %i",
                           build_gnomad, build_genome, an_aj))
}

miami <- don %>%
  plot_miami

snakemake@output[[1]] %>%
  ggsave(plot = miami, dpi = 600, width = 10, height = 7)
