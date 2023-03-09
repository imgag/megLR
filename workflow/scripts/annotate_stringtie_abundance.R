#!/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(annotables))

# Define default in/out files for testing if not run inside snakemake
if (exists("snakemake")) {
  abund <- snakemake@input$abund
  classif <- snakemake@input$classif
  out_tsv <- snakemake@output$tsv
  build <- snakemake@params$build
}

dta <- fread(abund)
dta <- dta[Coverage > 0,,]
dtc <- fread(classif)

dtc[,c("STRG_Gene_ID","STRG+Isoform_ID") := tstrsplit(isoform, "[.]", keep=c(2,3))]
dtc[,STRG_gene := paste0("STRG.",STRG_Gene_ID)]

dt <- merge(dta, dtc, by.x = "Gene ID", by.y =  "STRG_gene")
dt <- dt[,.(
  "chr" = Reference,
  "start" = Start,
  "end" = End,
  "strand" = Strand,
  "strg_isoform" = isoform,
  associated_gene,
  associated_transcript, 
  "coverage" = Coverage,
  TPM, FPKM, 
  length, exons, ref_length, ref_exons,
  structural_category,
  subcategory, 
  coding,
  ORF_seq
)]
dt[,ensid:=tstrsplit(associated_gene, '[.]', keep=1)]

dtm <- merge(
  dt, 
  as.data.table(grch38)[,.(ensgene, symbol, description),],
  by.x="ensid",
  by.y="ensgene",
  all.x=TRUE)
dtm <- dtm[,!"ensid",]
dtm <- dtm[order(chr, start)]

fwrite(dtm, out_tsv, sep="\t")
