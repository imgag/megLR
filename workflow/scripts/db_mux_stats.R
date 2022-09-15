library(ggplot2)
library(data.table)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# Define default in/out files for testing if not run inside snakemake
if (exists("snakemake")) {
  input_csv <- snakemake@input$csv
  out_mux <- snakemake@output$mux
} else {
  input_csv <- "test_qc_db/mux_scan_data_PAM61196_992887db.csv"
  out_mux <- "test_qc_db/out_mux.csv"
}

dt <- data.table::fread(input_csv)

s <- dt[,.(
  mux_scan_assessment = first(mux_scan_assessment), 
  total_good_samples = sum(total_good_samples), 
  total_time_active = sum(total_time_active),
  pore_median = sum(pore_median),
  pore_sd = sum(pore_sd),
  unavailable = sum(unavailable),
  multiple = sum(multiple),
  adapter = sum(adapter),
  unclassed = sum(unclassed),
  transition = sum(transition),
  zero = sum(zero)), 
  by=.(channel,scan_number,seconds_since_start_of_run)]


get_x <- function(i) {
  i <- i-1
  x <- ((i %/% 250) * 10 + (i %% 10)) +1
  return(x)
}

get_y <- function(i) {
  i <- i-1
  y <- ((i %/% 10) %% 25 ) +1
  return(y)
}

s[,x := get_x(channel),]
s[,y := get_y(channel),]

fwrite(s, out_mux)
