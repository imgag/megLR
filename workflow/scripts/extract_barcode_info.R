# Dual Demux analysis and copy

sink(file=snakemake@log[[1]], type = c("output"))
sink(file=NULL, type = c("message"))

print("Starting Dual Barcode Demux QC")

library(magrittr)

file <- snakemake@input[["demux_stats"]]
sample <- snakemake@wildcards[["sample"]]

# Read input file
print(paste("Reading input file", file))
diag <- readr::read_delim(file, delim = "\t", col_names = F ) %>%
  magrittr::set_names(c("sample","diag")) %>%
  dplyr::mutate(samples = purrr::map(stringr::str_split(stringr::str_trim(sample), " "), tail, -1)) %>%
  dplyr::mutate(bc_found = purrr::map(stringr::str_split(stringr::str_trim(diag), " "), head, 1))


# Print Sample Histogram
print("Counting sample occurences")
samplecount <- diag %>%
  dplyr::mutate(samples = ifelse (lengths(samples) > 1, "multiple", samples)) %>%
  tidyr::unchop(samples) %>%
  dplyr::count(samples) 

# Bar plot occurences
maxc <- samplecount %>% 
  dplyr::filter(samples != "unk" & samples != "multiple") %>%
  dplyr::summarise(max(n)) %>%
  dplyr::pull()

library(ggplot2)
sample_hist <- ggplot(samplecount, aes(x=samples, y=n)) +
  coord_flip(ylim=c(0,maxc)) +
  geom_bar(stat = "identity", color = "grey20", fill = "grey80") + 
  geom_text(aes(x="unk"),
            y = maxc,
            label = as.character(samplecount %>% dplyr::filter(samples == "unk") %>% dplyr::pull(n)),
            color = "grey10", hjust = "right", size = 3) +
  geom_text(aes(x="multiple"),
            y = maxc,
            label = as.character(samplecount %>% dplyr::filter(samples == "multiple") %>% dplyr::pull(n)),
            color = "grey10", hjust = "right", size = 3) +
  labs(y="Number of binned reads", x = "Samples") +
  theme_bw() 

ggsave(snakemake@output[["pdf"]], height = 5 + nrow(samplecount)*0.12)
ggsave(snakemake@output[["png"]], height = 5 + nrow(samplecount)*0.12)

# Export Files with multiple matches
diag %>%
  dplyr::filter(lengths(samples) > 1) %>%
  dplyr::pull(samples)


mult_files <- diag %>%
  dplyr::filter(lengths(samples) > 1) %>%
  dplyr::pull(samples) %>%
  stringi::stri_join_list(sep="_") %>%
  stringr::str_c("Sample_", sample, "/demux/",sample, "_", ., ".fastq")

mult_files <- c(mult_files, 
                paste0("Sample_", sample, "/demux/", sample, "_Multiple_Matches.fastq"),
                paste0("Sample_", sample, "/demux/", sample, "_unk.fastq"))

# Write file with multi sample FASTQ paths
write(mult_files, snakemake@output[["mult_files"]])

# Write file with correctly demuxed sample IDs
write(samplecount %>% dplyr::filter(samples != "unk" & samples != "multiple") %>% dplyr::pull(samples), snakemake@output[["demux_samples"]])