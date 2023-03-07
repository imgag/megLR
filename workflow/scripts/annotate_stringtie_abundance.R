#!/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(ggplot2))

# Define default in/out files for testing if not run inside snakemake
if (exists("snakemake")) {
  input_csv <- snakemake@input$csv
  out_device_pos <- snakemake@output$device_pos
  out_len_yield <- snakemake@output$len_yield
  out_timeline <- snakemake@output$timeline
} else {
  input_csv <- "table_out.csv"
  out_device_pos <- "plot_device_pos.png"
  out_len_yield <- "plot_len_yield.png"
  out_timeline <- "plot_timeline.png"
}
