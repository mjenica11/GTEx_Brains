# Write significant female- and male-biased genes in each region to file

# libraries
library(tidyverse)
library(mashr)

# constants
BASE <- "/scratch/mjpete11/human_monkey_brain/mashr/"

# Input
mash_results = readRDS(paste0(BASE, "output/mashr_results.rds"))
mash_beta = get_pm(mash_results)
mash_lfsr = get_lfsr(mash_results)

# Output
# csvs of sDEGs in each region
MALE <- "male_gene_counts/male_lst.rds"
FEMALE <- "female_gene_counts/female_lst.rds"

#_______________________________________________________________________________
# Make nested list to store genes that pass filtering threshold in each region
#_______________________________________________________________________________
# Significant if local false sign rate < 0.2 
lfsr_df <- rownames_to_column(as.data.frame(mash_lfsr), var = "genes")

# Move the gene column to the end
lfsr_df <- lfsr_df %>% relocate(genes, .after = last_col())

# Generate a nested list of genes that passed filter in each region
regions <- colnames(mash_lfsr)
keep_genes <- list()
for (i in seq(1:11)){
	keep_genes[[regions[[i]]]] <- with(lfsr_df, genes[lfsr_df[, i] < 0.2])
}

#_______________________________________________________________________________
# Split effect size matrix into two; positive and negative values
#_______________________________________________________________________________
# Replace values that do not meet subsetting condition with NA
male_upreg <- ifelse(mash_beta > 0, mash_beta, NA) 
female_upreg <- ifelse(mash_beta < 0, mash_beta, NA) 

# Convert df to list of lists
male_upreg <- apply(male_upreg, 2, as.list)
female_upreg <- apply(female_upreg, 2, as.list)

#_______________________________________________________________________________
# Subset sex-biased genes to inlcude only genes that passed the filtering
# threshold in either direction 
#_______________________________________________________________________________
male_filtered <- list()
for (i in seq(1:11)) {
	male_filtered[[i]] <- male_upreg[[i]][keep_genes[[i]]]
}

female_filtered <- list()
for (i in seq(1:11)) {
	female_filtered[[i]] <- female_upreg[[i]][keep_genes[[i]]]
}

# Drop genes that were significant but only in the other sex
male_filtered <- lapply(male_filtered, function(x) x[!is.na(x)])
female_filtered <- lapply(female_filtered, function(x) x[!is.na(x)])

# Name nested list by corresponding region
names(male_filtered) <- regions
names(female_filtered) <- regions

# Write nested lists to file
saveRDS(male_filtered, paste0(BASE, MALE))
saveRDS(female_filtered, paste0(BASE, FEMALE))

#_______________________________________________________________________________
# Reshape nested lists and write genes and betas for eacg region to a table
#_______________________________________________________________________________
# Make list of dfs from nested list
nested_lst_to_df <- function(nested_lst, strings) {
		dat <- stack(nested_lst, drop = FALSE)
		dat <- dat[, c(2,1)]
		col1 <- paste0(strings, "_genes")
		col2 <- paste0(strings, "_beta")
		colnames(dat) <- c(col1, col2)
		return(dat)
}

mlst <- Map(nested_lst_to_df, nested_lst = male_filtered, strings = regions)
flst <- Map(nested_lst_to_df, nested_lst = female_filtered, strings = regions)

# Write tables to file
sapply(names(mlst), function(x) write.csv(mlst[[x]], 
					file = paste0(BASE, "male_gene_counts/", x, ".csv")))
sapply(names(flst), function(x) write.csv(flst[[x]], 
					file = paste0(BASE, "female_gene_counts/", x, ".csv")))
