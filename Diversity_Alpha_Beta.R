## This code calculates the alpha-diversity and the beta-diversity from a phyloseq object.

## Load the required packages
library(ggplot2)
library(phyloseq)
library(microbiome)
library(knitr) 
library(dplyr)
library(reshape2)

## Load the data (Using publicly available data in this example)
## This dataset is part of the phyloseq library
data(GlobalPatterns)
GlobalPatterns    ## View the Phyloseq object

## 1. Alpha Diversity calculations
# Remove OTU's that are not present in any of the samples
dataset_pruned<- prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)
dataset_pruned ## View the Phyloseq object 

## Calculate the alpha diversity for indices of interest. Here calculating for Shannon and Inverse Simpson
alpha_diversity_table <- alpha(dataset_pruned, index = c("shannon", "diversity_inverse_simpson"))
kable(head(alpha_diversity_table)) ## View the top 6 rows 

## Extract the metadata from the phyloseq object
metadata <- meta(dataset_pruned)
kable(head(metadata))

## Add alpha diversity table to the metadata
metadata$Shannon <- alpha_diversity_table$diversity_shannon 
metadata$InverseSimpson <- alpha_diversity_table$diversity_inverse_simpson

## ggplot for Shannon diversity 
Fig2 <- ggplot(metadata, aes(x=SampleType, y=Shannon, fill=SampleType)) + 
  geom_point(position = "jitter") + 
  geom_boxplot(alpha=0.5)+ ylab("Shannon Diversity") +
  theme(legend.position="none",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), 
        axis.text.y=element_text(size=12)) +
  ggtitle ("Shannon Diversity by Sample Type")

Fig2

## Assuming the data for Shannon diversity is not normally distributed, 
## perform Kruskal-wallis test for Shannon Diversity by Sample Type 
Kruskal_test <- kruskal.test(Shannon ~ SampleType, data = metadata)
Kruskal_test ## check for p-value 

## If the p-value is significant, perform post-hoc Dunn test
library(dunn.test)
pairwise_comparisons <- dunn.test(metadata$Shannon, metadata$SampleType, method="bh",
                                  table = FALSE, kw = TRUE, list=TRUE)

## Save the values from the test as a dataframe
df_pairwise <- cbind.data.frame(pairwise_comparisons$Z, pairwise_comparisons$P, 
                                pairwise_comparisons$P.adjusted, pairwise_comparisons$comparisons)
colnames(df_pairwise) <- c("Z-Statistic", "P-Value", "Adjusted P-Value", "Group Comparison")

## 2. Beta Diversity calculations
## Pre-processing steps
## Remove taxa with less than 5 counts in a minimum of 10% of samples
taxa_filtered <- filter_taxa(GlobalPatterns, function(x) sum(x > 5) > (0.1*length(x)), TRUE)
taxa_filtered ## Number of taxa left after filtering

## Transform the sample counts for relative abundance 
dataset_transformed <- transform_sample_counts(taxa_filtered, function(x) x / sum(x))
# Preview the first 5 rows and columns of the transformed OTU table:
otu_table(dataset_transformed)[1:5, 1:5]

## Agglomerate taxa to Phylum level
dataset_Phylum <- tax_glom(dataset_transformed, taxrank = "Phylum")
ntaxa(dataset_Phylum) ## Number of taxa at Phylum level

## get the abundance table from the phyloseq object
abundance_table <- as.data.frame(dataset_Phylum@otu_table)

## calculate the bray-curtis distance matrix and transpose
bray_curtis_dist <- vegan::vegdist(t(abundance_table), method = "bray")

## Convert bray-curtis distance object into matrix format
bray_dist <- as.matrix(bray_curtis_dist) 
head(bray_dist)[,1:6] ##print the first 6 rows and columns

## Perform Dimensionality reduction using PCoA on bray-curtis distance matrix
pcoa <- ecodist::pco(bray_curtis_dist) 

## Store the values for first two principal coordinates and the sample_id's in a dataframe
pcoa_df <- data.frame(sample_id = row.names(pcoa$vectors),
                                  pcoa1 = pcoa$vectors[,1], 
                                  pcoa2 = pcoa$vectors[,2])

## Extract the eigenvalues for each principal coordinate 
## and calculate the percentage of total variance
eigenvals <- pcoa$values
var_explained <- as.vector(eigenvals/sum(eigenvals[eigenvals>0])*100) 

## Extract sample data from Phyloseq object
sample_data <- dataset_Phylum@sam_data

## Combine the sample data and the dataframe with principal coordinates for a ggplot
bray_curtis_pcoa_df <- cbind(pcoa_df,sample_data)
Fig3 <-ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = SampleType)) +
  geom_point() +
  labs(x=paste("PCoA1 (",round(var_explained[1],1),"%)"),
       y = paste("PCoA2 (",round(var_explained[2],1),"%)")) + theme_bw() + theme(legend.title=element_blank()) 

Fig3




