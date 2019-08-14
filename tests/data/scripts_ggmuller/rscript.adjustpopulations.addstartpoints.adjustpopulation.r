library("dplyr")

lag_gens <- function(x) {
	ans <- lag(all_gens_list)[which(all_gens_list == x)]
	if(is.na(ans)) return(0)
	return(ans)
}

adjustPopulation <- function(pop_df, first_gens, start_positions = 0.5) {
	
	# function to get the generation previous to a specified generation:

	
	# copy all rows for generations at which new genotypes appear:
	gens_list <- unique(first_gens$start_time)
	new_rows <- filter_(pop_df, ~Generation %in% gens_list)
	prev_rows <- filter_(pop_df, ~Generation %in% sapply(gens_list, lag_gens))
	# adjust generations of copied rows:
	new_rows$Generation <- new_rows$Generation - start_positions * (new_rows$Generation - sapply(new_rows$Generation, lag_gens))
	# adjust populations of copied rows:
	new_rows$Population <- (1 - start_positions) * new_rows$Population + start_positions * prev_rows$Population
	# add the copied rows to the dataframe:
	#pop_df <- bind_rows(pop_df, new_rows) %>%
	#	arrange_(~Generation, ~Identity)
	
	return(new_rows)
}

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
all_gens_list <- unique(pop_df$Generation)
first_gens = getInitialGenerations(pop_df)
result <- adjustPopulation(pop_df, first_gens)

write.table(result, args[2], sep = "\t", row.names = FALSE)