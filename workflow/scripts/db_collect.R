#!/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(openxlsx))

#print(str(snakemake))
folder <- snakemake@params$db_root

# Regex to match folder names
subf <- na.omit(str_match(list.dirs(paste0(folder, "/runs"), recursive = F), regex("(?:Q\\w{9}_)?\\d{5}\\w*_\\d{5}$")))

dt <- tibble(
  "start_time" = character(),
  "flow_cell_id" = character(),
  "project" = character(),
  "sample" = character(),
  "run_id" = character(),
  "qbic_id" = character(),
  "exp_type" = character(),
  "kit" = character(),
  "flow_cell_position" = character(),
  "sequencer" = character(),
  "flow_cell_product_code" = character(),
  "run_duration" = numeric(),
  "reads_number" = numeric(),
  "bases_number" = numeric(),
  "N50" = numeric())

for (f in subf) {

# Read report header
    report <- 
      read_lines(
      paste0(folder, "/runs/", f, "/", f, ".report.md"),
      skip = 4,
      n_max = 39,
      ) %>% 
      str_remove_all(regex("[\\\" ,]")) %>%
      str_split(":", n = 2, simplify = T) %>%
      as_tibble() 
  
  # Read pycoQC json
    qc_data <- rjson::fromJSON(file = paste0(folder, "/runs/", f, "/", f, ".pycoQC.json")) %>%
      purrr::pluck('Pass Reads') %>%
      purrr::flatten() %>%
      purrr::keep(names(.) %in% c(
          'run_duration',
          'reads_number',
          'bases_number',
          'N50'
          ))%>%
      tibble::as_tibble()
  
  # Extract data
  dt_rep <- report %>%
    spread(key = "V1", value = "V2") %>%
    select(c(
      "start_time" = "exp_start_time",
      "flow_cell_id",
      "project" = "protocol_group_id",
      "sample" = "sample_id",
      "exp_type" = "exp_script_name",
      "flow_cell_position" = "device_id",
      "sequencer" = "hostname",
      "flow_cell_product_code",
    ))
  
  # Extract info from folder name
  sf <- str_split(f, "_", simplify = T) 
  run_id <- sf[length(sf)]
  qbic_id <- if (length(sf) == 3) sf[1] else ""
  sample_id <- if (length(sf) == 3) sf[2] else sf[1]

  # Extract same info name from ONT report file
  s <- dt_rep %>% pull(sample)
  s <- str_split(s, "_", simplify = T) 
  run_report <- s[length(s)]
  qbic_id_report <- if (length(s) == 3) s[1] else ""
  sample_report <- if (length(s) == 3) s[2] else s[1]

  # Compare values and print mismatches
  if(run_id != run_report) {
    print(paste0('WARNING: run ID from folder (', run_id, ") does not match ID from ONT (", run_report, ')'))
  }
  if(qbic_id != qbic_id_report) {
    print(paste0('WARNING: QBIC ID from folder (', qbic_id, ") does not match ID from ONT (", qbic_id_report, ')'))
  }
  if(sample_id != sample_report) {
    print(paste0('WARNING: run sample name from folder (', sample_id, ") does not match ID from ONT (", sample_report, ')'))
  }

  dt_rep <- dt_rep %>%
    mutate(run_id = run_id) %>%
    mutate(sample = sample_id) %>%
    mutate(qbic_id = qbic_id)
  
  #  Split information from EXP type
  stype <- dt_rep %>% 
    pull(exp_type) %>%
    str_split("[:/]", simplify = T)

  s_exp <- if (length(stype) == 4) stype[[2]] else ""
  s_kit <- if (length(stype) == 4) stype[[4]] else ""
  
  dt_rep <- dt_rep %>%
    mutate(exp_type = s_exp) %>%
    mutate(kit = s_kit)

  # Add to main table
  dt <- add_row(dt, tibble_row(bind_cols(dt_rep, qc_data)))
}

print(dt)
# Export data to CSV
write_csv(dt, snakemake@output$csv)

dtt <- dt %>%
  mutate(start_time = lubridate::as_datetime(start_time)) %>%
  mutate(run_duration = lubridate::dhours(run_duration))

# Functions for Excel styling
wb <-
  openxlsx::createWorkbook()

sheet <-
  "Runs"
openxlsx::addWorksheet(wb, sheet)
openxlsx::writeDataTable(wb,
  sheet = sheet,
  x = dtt,
  tableStyle = "TableStyleMedium2",
  headerStyle = openxlsx::createStyle(
    halign = "center",
    valign = "center",
    wrapText = TRUE
  )
)
openxlsx::setColWidths(wb, sheet, cols = c(1:3, 5:6, 8:9), widths = 25)
openxlsx::saveWorkbook(wb, snakemake@output$xlsx, overwrite = TRUE)