# Make dendrograms based on total gene expression and sex biased genes

# libraries
library(mashr)
library(data.table)
library(dendextend)

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain/mashr/"

# Input
female <- readRDS(paste0(BASE, "female_gene_counts/female_lst.rds"))
male <- readRDS(paste0(BASE, "male_gene_counts/male_lst.rds"))
mash_results = readRDS(paste0(BASE, "output/mashr_results.rds"))

# Output
ALL_GENES <- paste0(BASE, "dendrograms/all_genes.pdf")
MALE <- paste0(BASE, "dendrograms/male.pdf")
FEMALE <- paste0(BASE, "dendrograms/female.pdf")

#_______________________________________________________________________________
# tissue phylogeny based on female- ande male-baised genes only  
#_______________________________________________________________________________
# Prep region names
regions <- c("amygdala", "anterior_cortex", "caudate", "cerebellum",
			 "frontal_cortex", "hippocampus", "hypothalamus", 
			 "nucleus_accumbens", "putamen", "spinal_cord", "substantia_nigra")

# female sDEG matrix
fem_df <- t(rbindlist(female, fill = TRUE))
colnames(fem_df) <- regions

# male sDEG matrix
male_df <- t(rbindlist(male, fill = TRUE))
colnames(male_df) <- regions

# Replace missing values with zero (for now)
fem_df[is.na(fem_df)] <- 0
male_df[is.na(male_df)] <- 0

# Similarity matrix
fem_dist <- dist(t(fem_df), method = "euclidean") 
male_dist <- dist(t(male_df), method = "euclidean") 

# Complete linkage clustering
fres <- hclust(fem_dist, method = "complete")
mres <- hclust(male_dist, method = "complete")

# Plot function
plot_dendrogram <- function(output, object, head_title) {
		pdf(output)
		par(cex = 0.8)
		plot(object, main = NULL, ylab = NULL, xlab = NULL)
		par(cex = 1)
		title(main = paste(head_title,  
	    "Euclidean distance with complete linkage clustering", sep = "\n"), 
	  	xlab = NULL, ylab = "Height")
		dev.off()
}
plot_dendrogram(output = FEMALE, object = fres, head_title = c("Dendrogram of significant female effect sizes;"))
plot_dendrogram(output = MALE, object = mres, head_title = c("Dendrogram of significant male effect sizes;"))

#_______________________________________________________________________________
# dendrogram of all gene effects 
#_______________________________________________________________________________
# Prepare df
gene_df <- get_pm(mash_results)
colnames(gene_df) <- regions

# Euclidean distance with complete, avergae, and single linkage clustering  
# Calculate the Euclidean distance across genes within the same region
genes <- dist(t(gene_df), method = "euclidean") 

# Complete linkage clustering
res1 <- hclust(genes, method = "complete")

# plot
plot_dendrogram(output = ALL_GENES, object = res1, head_title = c("Dendrogram of all effect sizes;"))
