#!/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(ggplot2))

#print(str(snakemake))
#print(snakemake)

# Load plot theme
snakemake@source("theme_cas.R")

# Load data
dt <- read_csv(snakemake@input$csv, col_types = cols())
dtt <- dt %>%
  mutate(start_time = lubridate::as_datetime(start_time)) %>%
  mutate(run_duration = lubridate::dhours(run_duration))

# _______________________
# Runs/Yield per week
dttw <- dtt %>%
  mutate(week = lubridate::floor_date(start_time, "week"))

p_weeks <- ggplot(dttw, aes(x = week, fill = project)) +
  geom_histogram(col = "grey20", binwidth = as.numeric(lubridate::dweeks(x=1))) +
  labs(x = "", y = "Flowcells", fill = "Project") +
  mytheme
ggsave(
  snakemake@output$timeline, 
  plot = p_weeks, 
  height = 5, 
  width = 8)

# ________________________
# Scatterplot length/yield
p_len_yield <- ggplot(dttw, aes(x = N50, y = bases_number, fill = project)) +
  geom_point(shape=21, col = "grey20", alpha = 0.9) + 
  mytheme +
  labs(
    x = "Read length (N50)",
    y = "Sequencing yield (bases)")
ggsave(
  snakemake@output$len_yield, 
  plot=p_len_yield, 
  height = 5, 
  width = 8)

# ________________________
# Per Device/Position
p_position <- 
  dttw %>% 
    filter(sequencer == "PC24B152") %>%
    separate(
      col = flow_cell_position,
      sep = 1,
      into = c("pos_row", "pos_col")) %>%
    group_by(pos_row, pos_col) %>%
    summarise(
      mean_yield=mean(bases_number),
      n_flowcells=n(),
      mean_length=mean(N50)) %>%
    gather(
      key = "key",
      value = "value",
      'mean_yield',
      'n_flowcells',
      'mean_length') %>%
    ggplot(
      aes(
        x = pos_row,
        y = pos_col,
        fill = value),
    ) + 
    geom_tile() +
    coord_fixed(ratio=1) +
    facet_grid(cols = vars(key)) + 
    mytheme 
ggsave(
  snakemake@output$device_pos, 
  plot = p_position, 
  height = 5, 
  width = 8)
