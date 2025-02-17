## This code calculates the relative abundance of Phylum level taxa from a phyloseq object and plots a bargraph using ggplot

## Load required libraries
library(phyloseq)
library(ggplot2)
library(RColorBrewer)

## Load the data (Using publicly available data in this example)
## This dataset is part of the phyloseq library
data(GlobalPatterns)

GlobalPatterns    ## View the Phyloseq object

## Pre-processing steps
# 1. Remove taxa with less than 5 counts in a minimum of 10% of samples
taxa_filtered <- filter_taxa(GlobalPatterns, function(x) sum(x > 5) > (0.1*length(x)), TRUE)
taxa_filtered ## Number of taxa left after filtering

# 2. Transform the sample counts for relative abundance 
dataset_transformed <- transform_sample_counts(taxa_filtered, function(x) x / sum(x))
# Preview the first 5 rows and columns of the transformed OTU table:
otu_table(dataset_transformed)[1:5, 1:5]

## Agglomerate taxa to Phylum level
dataset_Phylum <- tax_glom(dataset_transformed, taxrank = "Phylum")
ntaxa(dataset_Phylum) ## Number of taxa at Phylum level

## Select the top 10 Phylum level taxa
Top10 = names(sort(taxa_sums(dataset_Phylum), TRUE)[1:10])
bact_top_n_Phylum <- prune_taxa(Top10, dataset_Phylum) ## Phyloseq object with top 10 Phylum level taxa

## Convert to a data frame for ggplot using psmelt
rel_abundance_df <- psmelt(bact_top_n_Phylum)
rel_abundance_df$Phylum <- as.factor(rel_abundance_df$Phylum)

## Define customized color palette 
no_colors <- length(levels(rel_abundance_df$Phylum))
mycolors <- colorRampPalette(brewer.pal(11, "RdBu"))(no_colors)

Fig1 <-ggplot(rel_abundance_df, aes(x = SampleType, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bar chart
  scale_y_continuous(labels = scales::percent) +    
  labs(x = "Sample",
       y = "Relative Abundance") + scale_fill_manual(values = mycolors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Fig1

## save figure to local folder
ggsave("Fig1.pdf", plot = Fig1, width = 5, height = 4)




