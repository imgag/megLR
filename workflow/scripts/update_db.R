library(tidyverse)
library(tidyjson)
library(purrr)
library(rjson)

folder <- "/mnt/share/data/ont_run_qc"

subf <- na.omit(str_match(list.dirs(folder, recursive = F), regex("\\d{8}_\\d{4}_\\w+-\\w+-\\w+_\\D{3}\\d{5}_\\w{8}_\\w{5}$")))

for f in subf(
f <- subf[1]

# Read report header
report <- 
    read_lines(
    paste0(folder, "/", f, "/", f, ".report.md"),
    skip = 4,
    n_max = 39,
    ) %>% 
    str_remove_all(regex("[\\\" ,]")) %>%
    str_split(":", n = 2, simplify = T) %>%
    as_tibble()

# Read pycoQC json
rjson::fromJSON(file = paste0(folder, "/", f, "/", f, ".pycoQC.json")) %>%
    purrr::pluck('Pass Reads') %>%
    purrr::flatten() %>%
    purrr::keep(names(.) %in% c(
        'run_duration',
        'reads_number',
        'bases_number',
        'N50'
        ))%>%
    tibble::as_tibble() %>%
    print()
