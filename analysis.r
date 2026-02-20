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
genelevel_expr <- read.csv("data/GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv")
colnames(genelevel_expr)[1] <- c("ensembl_gene_id")

## Metadata
genelevel_metadata <- read.csv("data/GSE60450_filtered_metadata.csv")
colnames(genelevel_metadata)[c(1, 4)] <- c("sample_id", "developmental_stage")

# Inspect the data -------------------------

# What are the dimensions of each data set? (How many rows/columns in each?)
# Keep the code here for each file.


## Expression data
genelevel_expr %>% summarise(
  rows = n(),
  columns = ncol(.)
)

## Metadata
genelevel_metadata %>% summarise(
  rows = n(),
  columns = ncol(.)
)

# Prepare/combine the data for plotting ------------------------
# How can you combine this data into one data.frame?

## Transpose expression data by sample ID
genelevel_expr_long <- genelevel_expr %>%
  pivot_longer(cols = -c(ensembl_gene_id, gene_symbol),
               names_to = "sample_id",
               values_to = "expression")

## Check if sample IDs match
setequal(genelevel_expr_long$sample_id, genelevel_metadata$sample_id)

## Combine expression data and metadata by sample ID
genelevel_data <- genelevel_expr_long %>%
  left_join(genelevel_metadata, by = "sample_id")

# Plot the data --------------------------
## Plot the expression by cell type
## Can use boxplot() or geom_boxplot() in ggplot2

plot_genexpr <- function(data, gene_symbol) {
  
  # Filter for selected gene
  data_filtered <- data %>%
    filter(gene_symbol == !!gene_symbol)
  data_filtered$immunophenotype <- as.factor(data_filtered$immunophenotype)
  # Plot filtered data
  ggplot(data_filtered, 
         aes(x = immunophenotype, 
             y = expression,
             fill = immunophenotype)) +
    geom_boxplot(width = 0.4) +
    geom_jitter(position = position_jitter(width = 0.2),
                alpha = 0.6,
                size = 1) +
    labs(title = paste(gene_symbol, "expression by cell type"),
         x = "Cell type",
         y = "Expression (TMM-normalised CPM)") +
    theme_bw() +
    scale_fill_manual(values = c(
      "basal cell population" = "indianred2",
      "luminal cell population" = "deepskyblue2")) +
    theme(
      # x-axis
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(face = "bold"),
      # y-axis
      axis.text.y = element_text(size = 10, color = "black"),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      axis.title.y = element_text(face = "bold"),
      # Grid
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey92"),
      # Title
      title = element_text(face = "bold"),
      # Legend
      legend.position = "none"
    )
}

## Generate plot
plot1 <- plot_genexpr(genelevel_data, "Gnai3")

## Save the plot
### Show code for saving the plot with ggsave() or a similar function
ggsave(filename = "results/gene expression boxplot_gnai3.png",
       plot = plot1, 
       width = 2000,
       height = 2000,
       units = "px",
       dpi = 300)
