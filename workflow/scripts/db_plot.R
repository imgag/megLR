#!/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
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

# Setup plot theme
my_palette <- function(...) {  
    # Colour Brewer Set1
    colours <- c(
        "#E41A1C",
        "#377EB8",
        "#4DAF4A",
        "#984EA3",
        "#FF7F00",
        "#FFFF33",
        "#A65628",
        "#F781BF",
        "#999999"
        )
        
    function(n) {
        if (n == 0) {
            stop("Must request at least one colour from a hue palette.", 
                call. = FALSE)
        } else if (n<10) {
            sel_colours = colours[1:n]
        } else {
            sel_colours = colorRampPalette(colours)(n)
        }
        sel_colours
    }
    
}

scale_fill_custom <- function(..., aesthetics = "fill") {
    
    discrete_scale(aesthetics, "custom", my_palette(), ...)
}

# Load plot theme
mytheme <- function(){
    # Change main theme
    theme_bw() %+replace% 
    theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    ) 
}

# Load data
dt <- read_csv(input_csv, col_types = cols())
dtt <- dt %>%
  mutate(start_time = lubridate::as_datetime(start_time)) %>%
  mutate(run_duration = lubridate::dhours(run_duration)) %>%
  mutate(project = as.factor(project))

# _______________________
# Runs/Yield per week
dttw <- dtt %>%
  mutate(week = lubridate::floor_date(start_time, "week"))

p_weeks <- ggplot(dttw, aes(x = week, fill = project)) +
  geom_histogram(col = "grey20", binwidth = as.numeric(lubridate::dweeks(x=1))) +
  labs(x = "", y = "Flowcells", fill = "Project") +
#  scale_x_date(date_minor_breaks = "1 week") +
  mytheme() +
  scale_fill_custom()
ggsave(
  out_timeline,
  plot = p_weeks, 
  height = 5, 
  width = 8)

# ________________________
# Scatterplot length/yield
p_len_yield <- ggplot(dttw, aes(x = N50, y = bases_number, fill = project)) +
  geom_point(shape=21, col = "grey20", alpha = 0.9) + 
  mytheme() +
  scale_fill_custom() +
  labs(
    x = "Read length (N50)",
    y = "Sequencing yield (bases)")
ggsave(
  out_len_yield, 
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
    mytheme()
ggsave(
  out_device_pos, 
  plot = p_position, 
  height = 5, 
  width = 8)
