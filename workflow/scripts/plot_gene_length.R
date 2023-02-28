# Script to plot gene length distribution based on prokka prediction results

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
prokka_tsv <- args[1]
outdir <- args[2]
sample <- args[3]

# Read data
prokka_tsv_df <- read.csv(file=prokka_tsv, sep="\t", header=TRUE) %>%
  filter(ftype == "CDS")

prokka_tsv_df$length_bin <- addNA(cut(prokka_tsv_df$length_bp, 
                                  breaks=seq(0, 2000, 100)))
levels(prokka_tsv_df$length_bin) <- c("0-100", "100-200", "200-300", "300-400", "400-500",
                                      "500-600", "600-700", "700-800", "800-900", "900-1000",
                                      "1000-1100", "1100-1200", "1200-1300", "1300-1400", "1400-1500",
                                      "1500-1600", "1600-1700", "1700-1800", "1800-1900", "1900-2000",
                                      ">=2000")
count_df <- prokka_tsv_df %>%
  group_by(length_bin) %>%
  summarise(gene_count = n())

# Plot

color_count = length(unique(count_df$length_bin))
get_palette = colorRampPalette(brewer.pal(n = 12, "Paired"))

max_count = max(count_df$gene_count)
max_axis = pmax(100, 100 * ceiling(max_count / 100))

p <- ggplot(count_df, aes(x=length_bin, y=gene_count)) + 
  geom_bar(aes(fill = length_bin), stat = "identity", width= 0.8) +
  scale_fill_manual(values = get_palette(color_count)) + 
  geom_text(aes(label = gene_count), size = 3, hjust = 0.5, vjust = -0.2) +
  labs(x = "Gene length (bp)", 
       y = "Number of Genes", 
       title = paste(sample, ": Gene Length Distribution")) +
  scale_y_continuous(breaks = seq(0, max_axis + 50, 50), 
                     limits = c(0, max_axis + 50), 
                     expand = c(0,0)) +
  theme(panel.background = element_rect(fill="transparent", color = "transparent"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.5, color="lightgrey", linetype = "dotted"),
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.position = "none")

# Output plots
png_output = file.path(outdir, paste0(sample, "_gene_length.png"))
ggsave(filename = png_output, plot = p, width = 12, height = 8, dpi = 300, device = png)


