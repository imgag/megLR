suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# Define default in/out files for testing if not run inside snakemake
if (exists("snakemake")) {
  input_mux <- snakemake@input$csv
  plot_scan <- snakemake@output$plot_scan
  plot_total <- snakemake@output$plot_total
} else {
  # input_mux <- "test_qc_db/21103LRa019_02212.mux.csv"
  input_mux <- "test_qc_db/QIGMI027AM_21061LRa017L1_02192.mux.csv"
  plot_scan <- "test_qc_db/plot_scan.png"
  plot_total <- "test_qc_db/plot_total.png"
}

# Create empty plots for incomplete MUX files.
con <- file(input_mux) 
if (length(readLines(con)) < 100) {
    p_blank <- ggplot() + geom_blank()
    ggsave(plot_total, plot = p_blank, height = 5, width = 5)
    ggsave(plot_scan, plot = p_blank, height = 5, width = 5)
    quit(save="no")
}

s <- data.table::fread(input_mux, colClasses="")

ss <- s[,.(
  mux_scan_assessment = names(sort(table(mux_scan_assessment),decreasing=TRUE))[1], 
  total_good_samples = sum(as.numeric(total_good_samples)), 
  total_time_active = sum(as.numeric(total_time_active)),
  pore_median = sum(as.numeric(pore_median)),
  pore_sd = sum(as.numeric(pore_sd)),
  multiple = sum(as.numeric(multiple)),
  adapter = sum(as.numeric(adapter)),
  unclassed = sum(as.numeric(unclassed)),
  transition = sum(as.numeric(transition)),
  zero = sum(as.numeric(zero)),
  unavailable = sum(as.numeric(unavailable))),
  by=.(channel)]


get_x <- function(i) {
  i <- i-1
  x <- ((i %/% 250) * 10 + (i %% 10)) +1
  return(x)
}

get_y <- function(i) {
  i <- i - 1
  y <- ((i %/% 10) %% 25 ) +1
}  

ss[,x := get_x(as.numeric(channel)),]
ss[,y := get_y(as.numeric(channel)),]
s[,x := get_x(as.numeric(channel)),]
s[,y := get_y(as.numeric(channel)),]

blue = '#8fc5e6' #unavailable
green = '#00ff00' # green sequencing

# Plot overall productivity per pore
p_total <- ggplot(ss, aes(x= x, y = y)) +
  geom_tile(aes(fill = total_good_samples/1e6)) + 
  geom_vline(xintercept = c(30,60,90)) + 
  theme_classic() +
  coord_fixed(ratio=1) +
  scale_fill_gradient(low = "#dfe6eb", high = "#1db512" ) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "bottom")
              

# Plot pore status for all scan steps
p_scan <- ggplot(s, aes(x= x, y = y)) +
  geom_tile(aes(fill = mux_scan_assessment)) + 
  geom_vline(xintercept = c(30,60,90)) + 
  theme_classic() +
  facet_grid(rows = vars(seconds_since_start_of_run)) + 
  coord_fixed(ratio=1) + 
  theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         legend.position = "bottom")

# Write output files
ggsave(plot_total, plot = p_total, height = 7, width = 8)
ggsave(plot_scan, plot = p_scan, height = 80, width = 8, limitsize = FALSE)
