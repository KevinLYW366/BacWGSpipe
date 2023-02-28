#################################################################
# Script to plot COG annotation results                         #
#   count and plot the number of genes in each COG category     #
# by Yewei on 20221025                                          #
#################################################################

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
cog_annotation_file <- args[1]
cog_db_dir <- args[2]
sample <- args[3]
output_dir <- args[4]

# COG database file with descriptions of COG functional categories (fun-20.tab)
cog_fun_file = file.path(cog_db_dir, "fun-20.tab")
cog_fun_df <- read.csv(cog_fun_file, sep="\t", header=FALSE, row.names=NULL)
colnames(cog_fun_df) <- c("COG_functional_category", "color_code", "COG_functional_category_description")

# read COG annotation file
cog_annotation_df <- read.csv(cog_annotation_file, sep="\t", header=TRUE, row.names=NULL)
cog_annotation_plot_df <- cog_annotation_df %>%
  select(gene_id, COG_id, COG_functional_category) %>% 
  separate_rows(COG_functional_category, sep = "") %>%
  filter(COG_functional_category != "")
cog_annotation_plot_df <- merge(cog_annotation_plot_df, cog_fun_df, 
                                by = "COG_functional_category", all.x = TRUE) %>%
  group_by(COG_functional_category, COG_functional_category_description, color_code) %>%
  summarise(gene_count = n_distinct(gene_id)) %>%
  mutate(COG_functional_category_all = paste(COG_functional_category, 
                                             COG_functional_category_description,
                                             sep = ": ")) %>%
  mutate(color_code = paste0("#", color_code))

max_count = max(cog_annotation_plot_df$gene_count)
max_axis = pmax(100, 100 * floor(max_count / 100))

color_count = length(unique(cog_annotation_plot_df$COG_functional_category_all))
get_palette = colorRampPalette(brewer.pal(n = 12, "Paired"))

# plot
p1 <- ggplot(cog_annotation_plot_df, aes(x = COG_functional_category, y = gene_count)) +
  geom_bar(aes(fill = COG_functional_category_all), stat = "identity", width = 0.8) +
  scale_fill_manual(values = get_palette(color_count)) +
  geom_text(aes(label = gene_count), size = 3, hjust = 0.5, vjust = -0.2) +
  labs(x = "", 
       y = "Number of Genes", 
       title = paste0(sample, ": COG Annotation"),
       fill = "COG Functional Categories") +
  scale_y_continuous(breaks = seq(0, max_axis + 100, 50), 
                     limits = c(0, max_axis + 100), 
                     expand = c(0,0)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(panel.background = element_rect(fill="white", color = "transparent"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.5, color="lightgrey", linetype = "dotted"),
        axis.text.x = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"))

# output plots
png_output = file.path(output_dir, paste0(sample, "_cog.png"))
ggsave(filename = png_output, plot = p1, width = 12, height = 8, dpi = 300, device = png)
pdf_output = file.path(output_dir, paste0(sample, "_cog.pdf"))
ggsave(filename = pdf_output, plot = p1, width = 12, height = 8, dpi = 300, device = "pdf")
