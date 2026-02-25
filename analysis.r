# ---------------------------------------------------------

# Melbourne Bioinformatics Training Program

# This exercise to assess your familiarity with R and git. Please follow
# the instructions on the README page and link to your repo in your application.
# If you do not link to your repo, your application will be automatically denied.

# Leave all code you used in this R script with comments as appropriate.
# Let us know if you have any questions!


# You can use the resources available on our training website for help:
# Intro to R: https://mbite.org/intro-to-r
# Version Control with Git: https://mbite.org/intro-to-git/

# ----------------------------------------------------------

# Load libraries -------------------
# You may use base R or tidyverse for this exercise

# ex. library(tidyverse)

library(tidyverse)

# Load data here ----------------------
# Load each file with a meaningful variable name.

## Expression data
expr <- read.csv("data/GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv")
colnames(expr)[1] <- c("ensembl_gene_id") # matches with biomaRt notation

## Metadata
metadata <- read.csv("data/GSE60450_filtered_metadata.csv")
colnames(metadata)[c(1, 4)] <- c("sample_id", "developmental_stage")

# Inspect the data -------------------------

# What are the dimensions of each data set? (How many rows/columns in each?)
# Keep the code here for each file.


## Expression data
dim(expr) # rows = 23735, columns = 14

## Metadata
dim(metadata) # rows = 12, columns = 4

# Prepare/combine the data for plotting ------------------------
# How can you combine this data into one data.frame?

## Reshape dataframe
expr_long <- expr %>%
  pivot_longer(cols = -c(ensembl_gene_id, gene_symbol),
               names_to = "sample_id",
               values_to = "expression")

## Check if sample IDs match
setequal(expr_long$sample_id, metadata$sample_id)

## Match expression data and metadata by sample ID
expr_data <- expr_long %>%
  left_join(metadata, by = "sample_id")

# Plot the data --------------------------
## Plot the expression by cell type
## Can use boxplot() or geom_boxplot() in ggplot2

plot_genexpr <- function(data, gene_symbol) {
  
  # Filter for specified gene
  data_filtered <- data %>%
    filter(gene_symbol == !!gene_symbol)
  
  data_filtered$immunophenotype <- factor(data_filtered$immunophenotype,
    levels = c("basal cell population", "luminal cell population"),
    labels = c("basal cells", "luminal cells"))
  
  # Convert CPM to logCPM (log2 transformation)
  data_filtered$logCPM <- log2(data_filtered$expression + 1) # avoids plotting error with CPM values of 0
  
  # Boxplot
  ggplot(data_filtered, 
         aes(x = immunophenotype, 
             y = logCPM,
             fill = developmental_stage)) +
    geom_boxplot(width = 0.4,
                 position = position_dodge(width = 0.6)) +
    geom_jitter(aes(colour = developmental_stage),
      position = position_jitterdodge(jitter.width = 0, dodge.width = 0.6),
      alpha = 0.6,
      size = 1,
      show.legend = FALSE) +
    ## y-axis breaks
    scale_y_continuous(
      breaks = seq(floor(min(data_filtered$logCPM, na.rm = TRUE)),
      ceiling(max(data_filtered$logCPM, na.rm = TRUE)),
      by = 0.5)) +
    
    ## Axes and plot labels
    labs(title = bquote(italic(.(gene_symbol))),
         x = "Cell type",
         y = bquote("Normalised expression ("*log[2]*"CPM)"),
         fill = "Developmental stage") +
    
    ## Plot theme
    theme_bw() +
    scale_fill_manual(values = c(
      "18.5 day pregnancy" = "indianred2",
      "2 day lactation" = "deepskyblue2",
      "virgin" = "chartreuse3")) +
    scale_colour_manual(values = c(
      "18.5 day pregnancy" = "indianred2",
      "2 day lactation" = "deepskyblue2",
      "virgin" = "chartreuse3")) +
    theme(
      ### x-axis
      axis.text.x = element_text(size = 10, color = "black"),
      ### y-axis
      axis.text.y = element_text(size = 10, color = "black"),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      ### Grid
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey92"),
      ### Title
      plot.title = element_text(hjust = 0.5),
      ### Legend
      legend.position = "right")
}

## Generate plot
print(plot1 <- plot_genexpr(expr_data, "Gnai3"))
print(plot2 <- plot_genexpr(expr_data, "Tfe3"))

## Save the plot
### Show code for saving the plot with ggsave() or a similar function
ggsave(filename = "results/gene expression boxplot_gnai3.png",
       plot = plot1, 
       width = 2000,
       height = 2000,
       units = "px",
       dpi = 300)

ggsave(filename = "results/gene expression boxplot_Tfe3.png",
       plot = plot2, 
       width = 2000,
       height = 2000,
       units = "px",
       dpi = 300)
