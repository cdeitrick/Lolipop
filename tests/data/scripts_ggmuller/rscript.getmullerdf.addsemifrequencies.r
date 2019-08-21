library("dplyr")
addSemiFrequencies <- function(pop_df) {
	pop_df <- pop_df %>% group_by_(~Generation) %>% 
		mutate(Frequency = (Population / sum(Population)) / 2) %>%
		ungroup()
	pop_df$Population <- pop_df$Population / 2 # because of the duplication
	pop_df$Frequency[is.nan(pop_df$Frequency)] <- 0
	
	return(pop_df)
}

args = commandArgs(trailingOnly=TRUE)
filename <- args[1]
pop_df <- read.table(filename, sep = "\t", header = TRUE)

result <- addSemiFrequencies(pop_df)

write.table(result, args[2], sep = "\t")