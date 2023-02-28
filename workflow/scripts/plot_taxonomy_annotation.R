#########################################################
# Script to plot eggNOG Taxonomy annotation results     #
#   count and plot the number of genes in each taxonomy #
#   based on NCBI Taxonomy database                     #
# by Yewei on 20221026                                  #
#########################################################

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
eggnog_tax_annotation_file <- args[1]
sample <- args[2]
output_dir <- args[3]

# read eggnog annotation file
eggnog_tax_annotation_df <- read.csv(eggnog_tax_annotation_file, sep="\t", header=TRUE, row.names=NULL)
eggnog_tax_plot_df <- eggnog_tax_annotation_df %>%
  group_by(NCBI_Taxonomy_id, NCBI_Taxonomy_name) %>%
  summarise(gene_count = n_distinct(gene_id)) %>%
  arrange(desc(gene_count)) %>%
  head(20) %>%
  mutate(NCBI_Taxonomy_id_name = paste(NCBI_Taxonomy_id, NCBI_Taxonomy_name, sep=": "))
levels = eggnog_tax_plot_df$NCBI_Taxonomy_id_name
eggnog_tax_plot_df$NCBI_Taxonomy_id_name = factor(eggnog_tax_plot_df$NCBI_Taxonomy_id_name, 
                                                  levels = levels)
levels = eggnog_tax_plot_df$NCBI_Taxonomy_id
eggnog_tax_plot_df$NCBI_Taxonomy_id = factor(eggnog_tax_plot_df$NCBI_Taxonomy_id, 
                                             levels = levels)

max_count = max(eggnog_tax_plot_df$gene_count)
max_axis = pmax(100, 100 * ceiling(max_count / 100))

color_count = length(unique(eggnog_tax_plot_df$NCBI_Taxonomy_id_name))
get_palette = colorRampPalette(brewer.pal(n = 12, "Paired"))

# plot
p1 <- ggplot(eggnog_tax_plot_df, aes(x = NCBI_Taxonomy_id, y = gene_count)) +
  geom_bar(aes(fill = NCBI_Taxonomy_id_name), stat = "identity", width = 0.8) +
  scale_fill_manual(values = get_palette(color_count)) +
  geom_text(aes(label = gene_count), size = 3, hjust = 0.5, vjust = -0.2) +
  labs(x = "Top20 NCBI Taxonomy IDs", 
       y = "Number of Genes", 
       title = paste0(sample, ": NCBI Taxonomy Annotation"),
       fill = "NCBI Taxonomy Names") +
  scale_y_continuous(breaks = seq(0, max_axis + 100, 500), 
                     limits = c(0, max_axis + 100), 
                     expand = c(0,0)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color = "transparent"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.5, color="lightgrey", linetype = "dotted"),
        axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"))

# output plots
png_output = file.path(output_dir, paste0(sample, "_taxonomy.png"))
ggsave(filename = png_output, plot = p1, width = 12, height = 8, dpi = 300, device = png)
pdf_output = file.path(output_dir, paste0(sample, "_taxonomy.pdf"))
ggsave(filename = pdf_output, plot = p1, width = 12, height = 8, dpi = 300, device = "pdf")
