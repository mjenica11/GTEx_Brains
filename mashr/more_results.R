#!/usr/bin/Rscript

# Purpose: Summarize mashr results
# libraries
library(tidyverse)
library(mashr)

# constants
BASE <- "/scratch/mjpete11/human_monkey_brain/mashr/"

# Input
female <- readRDS(paste0(BASE, "female_gene_counts/female_lst.rds"))
male <- readRDS(paste0(BASE, "male_gene_counts/male_lst.rds"))

#mash_results = readRDS(paste0(BASE, "output/mashr_results.rds"))
#mash_beta = get_pm(mash_results)
#mash_lfsr = get_lfsr(mash_results)
#
# Output
FMAJ <- paste0(BASE, "new_output/female_majority.csv")
MMAJ <- paste0(BASE, "new_output/male_majority.csv")
FMIN <- paste0(BASE, "new_output/female_unique.csv")
MMIN <- paste0(BASE, "new_output/male_unique.csv")
MHIST <- paste0(BASE, "new_output/male_histogram.pdf")
FHIST <- paste0(BASE, "new_output/female_histogram.pdf")
ALL_HIST <- paste0(BASE, "new_output/all_significant_histogram.pdf")

#_____________________________________________________________________________ 
# Which genes are sex-biased across a majority of regions (>=8)? 
#_____________________________________________________________________________ 
# Generate list of just the gene names in each region
count_genes <- function(lst_parameter, number) {
		llst <- vector("list", length(lst_parameter))
		nams <- names(lst_parameter)
		for (i in seq_along(lst_parameter)) {
				llst[[i]] <- names(lst_parameter[[i]])
		}
		res <- as.data.frame(table(unique(stack(setNames(llst, seq_along(llst))))$values))
		if (number == 1) {
				freq <- res[res$Freq == number, ]
		}
		else if (number > 1 | number == 0) {
				freq <- res[res$Freq >= number, ]
		}
		else { }
		# Rename cols and add row number
		colnames(freq) <- c("gene", "frequency")
		row.names(freq) <- seq(1, nrow(freq), 1)
		return(freq)
}
# Frequency of genes that are sig in >= 8 regions
mmaj <- count_genes(lst_parameter = male, number = 8)
fmaj <- count_genes(lst_parameter = female, number = 8)

# " " " " " " " any number of regions
mmin <- count_genes(lst_parameter = male, number = 1)
fmin <- count_genes(lst_parameter = female, number = 1)

write.csv(fmaj, FMAJ)
write.csv(mmaj, MMAJ)
write.csv(fmin, FMIN)
write.csv(mmin, MMIN)

#_____________________________________________________________________________ 
# Are sex-biased genes more likely to be region-specific or shared? 
#_____________________________________________________________________________ 
# Get frequency of all significant gene
mall <- count_genes(lst_parameter = male, number = 0)
fall <- count_genes(lst_parameter = female, number = 0)

# frequency of all sex-biased genes
all_genes <- rbind(mall, fall) 

# Histogram of the number of sex-biased genes vs number of regions 
# plot function
hist_function <- function(adress, object, sex) {
	pdf(adress)
	hist(object$frequency, 
		main = paste0("Distribution of significant ", sex, "-biased genes across n regions"),
	 	xlab = "Number of regions")
	dev.off()
}
hist_function(adress = MHIST, object = mall, sex = "male")
hist_function(adress = FHIST, object = fall, sex = "female")
hist_function(adress = ALL_HIST, object = all_genes, sex = "sex")
