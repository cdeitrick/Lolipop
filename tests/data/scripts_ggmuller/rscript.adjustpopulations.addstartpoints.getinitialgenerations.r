library("dplyr")

getInitialGenerations <- function(pop_df) {
	# set small initial population size:
	init_size <- 0
	
	# get reference list of generations at which new genotypes appear (and previous generations):
	min_gen <- min(pop_df$Generation)
	first_gens <- group_by_(pop_df, ~Identity) %>%
		filter_(~max(Population) > 0) %>% 
		summarise_(start_time = ~min(Generation[which(Population > 0)]), 
				   previous_time = ~lag(Generation)[min(which(Population > 0))]) %>%
		filter_(~start_time > min_gen) %>%
		ungroup()
	return(first_gens)
}
	
args = commandArgs(trailingOnly=TRUE)
#filename <- "/home/cld100/Documents/github/muller_diagrams/tests/test_get_muller_df_migration/data_add_start_points/B1_Muller.ggmuller.populations.tsv"
filename <- args[1]
pop_df <- read.csv(filename, sep = "\t", header = TRUE)
result <- getInitialGenerations(pop_df)

write.table(result, args[2], sep = "\t")
