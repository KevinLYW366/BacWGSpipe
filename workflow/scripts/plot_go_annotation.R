#################################################################
# Script to plot GO annotation results                          #
#   count and plot the number of genes in each GO level-1 class #
# by Yewei on 20221024                                          #
#################################################################

# Load packages
library(ggplot2)
library(dplyr)
library(ggtext)
library(glue)

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
go_annotation_file <- args[1]
sample <- args[2]
output_dir <- args[3]

# Read GO annotation results file
go_annotation_df <- read.csv(go_annotation_file, header = TRUE, sep = "\t", row.names = NULL) %>%
  filter(GO_level == 1) %>%
  group_by(GO_class, GO_description) %>%
  summarise(gene_count = n_distinct(gene_id)) %>%
  mutate(GO_description_color = 
           ifelse(GO_class == "biological_process", glue("<i style='color:#E41A1C'>{GO_description}</i>"), 
                  ifelse(GO_class == "cellular_component", glue("<i style='color:#377EB8'>{GO_description}</i>"),
                         glue("<i style='color:#4DAF4A'>{GO_description}</i>"))))

max_count = max(go_annotation_df$gene_count)
max_axis = pmax(100, 100 * ceiling(max_count / 100))

color_go <- c("biological_process" = "#E41A1C",
              "cellular_component" = "#377EB8",
              "molecular_function" = "#4DAF4A")

# plot
p1 <- ggplot(go_annotation_df, aes(x = GO_description_color, y = gene_count)) +
  geom_bar(aes(fill = GO_class), stat = "identity", width = 0.8) +
  geom_text(aes(label = gene_count), size = 2.3, hjust = 0.5, vjust = -0.3) +
  scale_fill_manual(values = color_go) +
  facet_grid(~ GO_class, 
             scales = "free", 
             space = "free") +
  labs(x = "", 
       y = "Number of Genes", 
       title = paste0(sample, ": GO Annotation (Level-2)"),
       fill = "GO Class") +
  scale_y_continuous(breaks = seq(0, max_axis + 100, 100),
                     limits = c(0, max_axis + 100), 
                     expand = c(0,0)) + 
  theme(panel.background = element_rect(fill="white", color = "transparent"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.5, color="lightgrey", linetype = "dotted"),
        panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_markdown(face = "bold", angle = 60, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.text.y = element_text(face = "bold"), 
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        strip.text.x = element_blank())

# output plots
png_output <- file.path(output_dir, paste0(sample, "_go.png"))
ggsave(filename = png_output, plot = p1, width = 13, height = 8, dpi = 300, device = png)
pdf_output <- file.path(output_dir, paste0(sample, "_go.pdf"))
ggsave(filename = pdf_output, plot = p1, width = 13, height = 8, dpi = 300, device = "pdf")

