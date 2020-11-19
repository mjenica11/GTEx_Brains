# Write significant female- and male-biased genes in each region to file

# libraries
library(tidyverse)
library(mashr)

# Input
mash_results = readRDS('/scratch/mjpete11/human_monkey_brain/mashr/output/mashr_results.rds')
mash_beta = get_pm(mash_results)
mash_lfsr = get_lfsr(mash_results)

# Output
FEMALE_UPREG <- "/scratch/mjpete11/human_monkey_brain/mashr/output/female_upregulated.csv"
MALE_UPREG <- "/scratch/mjpete11/human_monkey_brain/mashr/output/male_upregulated.csv"

#_______________________________________________________________________________
# Generate tables listing the female- and male-biased genes
#_______________________________________________________________________________
# Significant if local false sign rate < 0.2 
# Make nested list to store genes that pass filtering threshold in each region
lfsr_df <- rownames_to_column(as.data.frame(mash_lfsr), var = "rowname")
lfsr_df <- lfsr_df %>% relocate(rowname, .after = last_col())
regions <- colnames(mash_lfsr)
keep_genes <- list()
for (i in seq(1:11)){
	keep_genes[[regions[[i]]]] <- with(lfsr_df, rowname[lfsr_df[, i] < 0.2])
}

# Split effect size matrix into two; positive and negative values
# Replace values that do not meet subsetting condition with NA
male_upreg <- ifelse(mash_beta > 0, mash_beta, NA) 
female_upreg <- ifelse(mash_beta < 0, mash_beta, NA) 

# Convert df to list of lists
male_upreg <- apply(male_upreg, 2, as.list)
female_upreg <- apply(female_upreg, 2, as.list)

# Subset sex-biased genes to inlcude only genes that passed the filtering threshold
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

# Convert nested list to df
male_filtered <- sapply(male_filtered, '[', seq(max(sapply(male_filtered,length))))
female_filtered <- sapply(female_filtered, '[', seq(max(sapply(female_filtered,length))))

# Write to file
write.csv(male_filtered, MALE_UPREG)
write.csv(female_filtered, FEMALE_UPREG)
