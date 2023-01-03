#!/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))

# Define default in/out files for testing if not run inside snakemake
if (exists("snakemake")) {
  input_vcf <- snakemake@input$vcf
  out_device_pos <- snakemake@output$device_pos
  out_len_yield <- snakemake@output$len_yield
  out_timeline <- snakemake@output$timeline
} else {
  input_vcf <- "../../resources/testdata/repeat_expansion/example.repeats.vcf.gz"
  out_device_pos <- "plot_device_pos.png"
  out_len_yield <- "plot_len_yield.png"
  out_timeline <- "plot_timeline.png"
  setwd("~/Documents/Arbeit/IMGAG/dev/ont_tools/workflow/scripts")
}

vcf <- read_tsv(input_vcf, comment="##")

vcf %>% 
  pull(INFO) %>%
  str_split("[;/]", simplify = T) 

%>%
  str_split("[=]", simplify = T)
  