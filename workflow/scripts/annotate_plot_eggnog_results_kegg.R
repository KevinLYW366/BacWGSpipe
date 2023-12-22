####################################################################
# Script to annotate and plot KEGG annotation results              #
#   count and plot the number of genes in each KEGG level2 pathway #
# by Yewei on 20221018                                             #
####################################################################

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(gtable)
library(grid)
library(RColorBrewer)
library(ggtext)
library(glue)

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
eggnog_results_file <- args[1]
kegg_path_db_file <- args[2]
sample <- args[3]
output_dir <- args[4]


# convert eggnog results to a list of kegg pathway IDs (mapXXXXX)
eggnog_results_df <- read.csv(eggnog_results_file, comment.char = "#", header = FALSE, sep = "\t")
eggnog_kegg_path <- eggnog_results_df[13]
colnames(eggnog_kegg_path) <- c("pathway_id")
eggnog_kegg_path_list <- eggnog_kegg_path %>%
  separate_rows(pathway_id, sep = ",") %>%
  filter(pathway_id != "-" & !grepl('ko', pathway_id))

# read kegg reference pathway database
kegg_path_db_df <- read.csv(kegg_path_db_file, header = TRUE, sep = "\t")
colnames(kegg_path_db_df)[1] <- "pathway_id"

# merge two dataframes
# text output
df_tmp <- eggnog_results_df[c(1,12,13)]
colnames(df_tmp) <- c("gene_id", "KO_id", "pathway_id")

df_tmp <- df_tmp %>%
  separate_rows(pathway_id, sep = ",") %>% 
  filter(!grepl('ko', pathway_id)) %>%
  filter(KO_id != "-")

kegg_path_df_output <- merge(df_tmp, kegg_path_db_df, by = "pathway_id", all.x = TRUE) %>%
  select(-identified_KOs) %>%
  select(gene_id, KO_id, pathway_id, pathway_level1, pathway_level2, pathway_level3) %>%
  arrange(gene_id, pathway_id, pathway_level1, pathway_level2, pathway_level3)

kegg_path_df_output$pathway_level1[kegg_path_df_output$pathway_id == "-"] <- "-"
kegg_path_df_output$pathway_level2[kegg_path_df_output$pathway_id == "-"] <- "-"
kegg_path_df_output$pathway_level3[kegg_path_df_output$pathway_id == "-"] <- "-"

kegg_path_df_output <- kegg_path_df_output %>%
  filter(!is.na(pathway_level1)) %>%
  mutate(KO_id = str_replace_all(KO_id, "ko:", ""))

write.table(kegg_path_df_output, 
            file = paste0(output_dir, "/", sample, "_kegg.xls"),
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

# plot output
kegg_path_df_plot <- merge(eggnog_kegg_path_list, kegg_path_db_df, by = "pathway_id", all.x = TRUE)

# plot_data
plot_data <- kegg_path_df_plot %>%
  filter(!is.na(pathway_level1)) %>%
  group_by(pathway_level1, pathway_level2) %>%
  summarise(gene_count = n()) %>%
  arrange(pathway_level1, pathway_level2)

# color palette Set1
plot_data$color <- ifelse(plot_data$pathway_level1 == "Metabolism", "#FF7F00",
                          ifelse(plot_data$pathway_level1 == "Genetic Information Processing", "#4DAF4A", 
                                 ifelse(plot_data$pathway_level1 == "Environmental Information Processing", "#377EB8", 
                                        ifelse(plot_data$pathway_level1 == "Cellular Processes", "#E41A1C", 
                                               ifelse(plot_data$pathway_level1 == "Organismal Systems", "#A65628", 
                                                      ifelse(plot_data$pathway_level1 == "Human Diseases", "#984EA3", "#FFFF33")))))) 

# modify plot data
plot_data <- plot_data %>%
  mutate(
    pathway_level2_color = glue("<i style='color:{color}'>{pathway_level2}</i>"), 
    pathway_level1_color = glue("<i style='color:{color}'>{pathway_level1}</i>"))

color_pathway_level1 <- c("Cellular Processes" = "#E41A1C",
           "Environmental Information Processing" = "#377EB8", 
           "Genetic Information Processing" = "#4DAF4A", 
           "Human Diseases" = "#984EA3", 
           "Metabolism" = "#FF7F00", 
           "Organismal Systems" = "#A65628",
           "Drug Development" = "#FFFF33")

max_count = max(plot_data$gene_count)
max_axis = pmax(100, 100 * round(max_count / 100))

# plot
p1 <- ggplot(plot_data, aes(x = gene_count, y = pathway_level2_color)) +
  geom_bar(aes(fill = pathway_level1), stat = "identity") +
  scale_fill_manual(values = color_pathway_level1) +
  geom_text(aes(label = gene_count), size = 3, hjust = 0, vjust = 0.05) + 
  facet_grid(pathway_level1 ~ ., 
             scales = "free_y", 
             space = "free_y") +
  labs(x = "Number of Genes", 
       y = "", 
       title = paste0(sample, ": KEGG Pathway Annotation")) +
  scale_x_continuous(breaks = seq(0, max_axis, 100), 
                     limits = c(0, max_axis + 100), 
                     expand = c(0,0)) +
  theme(panel.spacing = unit(1, 'lines'),
        strip.text.y = element_markdown(angle = 0, size = 11, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", color = "transparent"), 
        legend.position="none", 
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_markdown(face = "bold"),
        axis.text.x = element_text(face = "bold", angle = 60, hjust = 1, vjust = 1),
        axis.title.x = element_text(face = "bold"),
        axis.line.x.bottom = element_line(color = "black"), 
        panel.background = element_rect(fill="white", color = "transparent"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.5, color="lightgrey", linetype = "dotted"))

# following codes set "space = 'free'" in facet_wrap like facet_grid
gt <- ggplotGrob(p1)
panels <-c(subset(gt$layout, grepl("panel", gt$layout$name), se=t:r))
for(i in rev(panels$t-1)) {
  gt = gtable_add_rows(gt, unit(0.5, "lines"), i)
}
panels <-c(subset(gt$layout, grepl("panel", gt$layout$name), se=t:r))
strips <- c(subset(gt$layout, grepl("strip-r", gt$layout$name), se=t:r))
stripText = gtable_filter(gt, "strip-r")
for(i in 1:length(strips$t)) {
  gt = gtable_add_grob(gt, stripText$grobs[[i]]$grobs[[1]], t=panels$t[i]-1, l=5)
}
gt = gt[,-6]
for(i in panels$t) {
  gt$heights[i-1] = unit(0.8, "lines")
  gt$heights[i-2] = unit(0.2, "lines")
}

# output plots
png_output = file.path(output_dir, paste0(sample, "_kegg.png"))
ggsave(filename = png_output, plot = gt, width = 10, height = 8, dpi = 100, device = png)
pdf_output = file.path(output_dir, paste0(sample, "_kegg.pdf"))
ggsave(filename = pdf_output, plot = gt, width = 10, height = 8, dpi = 300, device = "pdf")

unlink("Rplots.pdf")