# Make dendrograms based on total gene expression and sex biased genes

# libraries
library(tidyverse)
library(mashr)
library(stringr)

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain/mashr/"

# Input
mash_results = readRDS('/scratch/mjpete11/human_monkey_brain/mashr/output/mashr_results.rds')
mash_beta = get_pm(mash_results)
mash_sbet = get_pm(mash_results) / get_psd(mash_results)
mash_lfsr = get_lfsr(mash_results)

# Output
#ALL_EUCLIDEAN_COMPLETE <- paste0(BASE, "dendrograms/", "all_euclidean_complete.pdf")

#_______________________________________________________________________________
# Read in an organize sex-biased gene effect sizes in each region  
#_______________________________________________________________________________
regions <- colnames(mash_beta)

# Read in male results as list of dfs
setwd(paste0(BASE, "male_gene_counts/"))
mfiles <- list.files()
mlst <- lapply(mfiles, read.csv)

# Read in female " " " "
setwd(paste0(BASE, "female_gene_counts/"))
ffiles <- list.files()
flst <- lapply(ffiles, read.csv)

# Drop rownumber column and drop "." in column names
drop_col <- function(x) {
		x <- x[, -c(1)]
		colnames(x) <- str_replace_all(colnames(x), pattern = "\\.", "") 
		return(x)
}
mlst <- lapply(mlst, drop_col)
flst <- lapply(flst, drop_col)

#_______________________________________________________________________________
# Experiment with different similarty measures and clustering methods 
#_______________________________________________________________________________
# Euclidean distance with complete, avergae, and single linkage clustering  
# Calculate the Euclidean distance across genes within the same region
euclidean <- dist(t(mash_beta), method = "euclidean") 

# Complete linkage clustering
res1 <- hclust(euclidean, method = "complete")

# Plot cluster dendrogram
pdf(paste0(BASE, "dendrograms/", "all_euclidean_complete.pdf"))
plot(res1, main = "Euclidean distance with complete linkage clustering")
par(crt = 0.7)
dev.off()

# Euclidean ditance with average linkage clustering
res2 <- hclust(euclidean, method = "average")

# Plot cluster dendrogram
pdf(paste0(BASE, "dendrograms/", "all_euclidean_average.pdf"))
plot(res2, main = "Euclidean distance with average linkage clustering")
dev.off()

# Euclidean distance with single linkage clustering
res3 <- hclust(euclidean, method = "single")

# Plot cluster dendrogram
pdf(paste0(BASE, "dendrograms/", "all_euclidean_single.pdf"))
plot(res3, main = "Euclidean distance with average linkage clustering")
dev.off()

#_______________________________________________________________________________
# Manhattan distance with complete, avergae, and single linkage clustering 
#_______________________________________________________________________________
# Calculate the Euclidean distance between gene effect sizes
manhattan <- dist(t(mash_beta), method = "manhattan") 

# Complete linkage clustering
res4 <- hclust(manhattan, method = "complete")

# Plot cluster dendrogram
pdf(paste0(BASE, "dendrograms/", "all_manhattan_complete.pdf"))
plot(res4, main = "Manhattan distance with complete linkage clustering")
dev.off()

# Manhattan distance with average linkage clusterin
res5 <- hclust(manhattan, method = "average")

# Plot cluster dendrogram
pdf(paste0(BASE, "dendrograms/", "all_manhattan_average.pdf"))
plot(res5, main = "Manhattan distance with average linkage clustering")
dev.off()

# Manhattan distance with single linkage clustering
res6 <- hclust(manhattan, method = "single")

# Plot cluster dendrogram
pdf(paste0(BASE, "dendrograms/", "all_manhattan_single.pdf"))
plot(res6, main = "Manhattan distance with single linkage clustering")
dev.off()

#_______________________________________________________________________________
# Repeat with standardized betas 
#_______________________________________________________________________________
euclidean <- dist(t(mash_sbet), method = "euclidean") 

# Complete linkage clustering
res1 <- hclust(euclidean, method = "complete")

# Plot cluster dendrogram
pdf(paste0(BASE, "dendrograms/", "sbet_euclidean_complete.pdf"))
plot(res1, main = "Euclidean distance of standardized betas with complete linkage clustering")
dev.off()

# Euclidean ditance with average linkage clustering
res2 <- hclust(euclidean, method = "average")

# Plot cluster dendrogram
pdf(paste0(BASE, "dendrograms/", "sbet_euclidean_average.pdf"))
plot(res2, main = "Euclidean distance of standardized betas with average linkage clustering")
dev.off()

# Euclidean distance with single linkage clustering
res3 <- hclust(euclidean, method = "single")

# Plot cluster dendrogram
pdf(paste0(BASE, "dendrograms/", "sbet_euclidean_single.pdf"))
plot(res3, main = "Euclidean distance of standardized betas with average linkage clustering")
dev.off()
