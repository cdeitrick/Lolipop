library("dplyr")

add_start_points <- function(pop_df, start_positions = 0.5) {
	
	original_colname <- "Generation"
	# rename Time column (original name will be restored later):
	if("Time" %in% colnames(pop_df) && !("Generation" %in% colnames(pop_df))) {
		colnames(pop_df)[colnames(pop_df) == "Time"] <- "Generation"
		original_colname <- "Time"
	}
	
	# set small time interval:
	all_gens_list <- unique(pop_df$Generation)
	delta <- abs(min(1E-2 * min(diff(all_gens_list)), 1E-4 * (max(all_gens_list) - min(all_gens_list))))
	delta_debug <- min(diff(all_gens_list))
	start_positions <- max(start_positions, delta)
	start_positions <- min(start_positions, 1 - delta)
	start_positions <- max(start_positions, 0.5)
	
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
	
	# if all genotypes appear at the first time point then don't make any changes:
	if(dim(first_gens)[1] == 0) return(pop_df)
	
	# function to get the generation previous to a specified generation:
	lag_gens <- function(x) {
		ans <- lag(all_gens_list)[which(all_gens_list == x)]
		if(is.na(ans)) return(0)
		return(ans)
	}
	
	# copy all rows for generations at which new genotypes appear:
	gens_list <- unique(first_gens$start_time)
	new_rows <- filter_(pop_df, ~Generation %in% gens_list)
	prev_rows <- filter_(pop_df, ~Generation %in% sapply(gens_list, lag_gens))
	# adjust generations of copied rows:
	new_rows$Generation <- new_rows$Generation - start_positions * (new_rows$Generation - sapply(new_rows$Generation, lag_gens))
	# adjust populations of copied rows:
	new_rows$Population <- (1 - start_positions) * new_rows$Population + start_positions * prev_rows$Population
	# add the copied rows to the dataframe:
	pop_df <- bind_rows(pop_df, new_rows) %>%
		arrange_(~Generation, ~Identity)
	
	# adjust generations in reference list:
	first_gens$Generation <- first_gens$start_time - start_positions * (first_gens$start_time - sapply(first_gens$start_time, lag_gens))
	# set small initial populations in reference list:
	first_gens$Population2 <- init_size
	
	# replace initial populations in the dataframe with values from reference list:
	pop_df <- merge(pop_df, first_gens, all.x = TRUE)
	pop_df$Population <- ifelse(is.na(pop_df$Population2), pop_df$Population, pop_df$Population2)
	pop_df <- pop_df[, !(names(pop_df) %in% c("Population2", "start_time", "previous_time"))]
	
	# restore original time column name:
	colnames(pop_df)[colnames(pop_df) == "Generation"] <- original_colname
	
	return(pop_df)
}
args = commandArgs(trailingOnly=TRUE)
#filename <- "/home/cld100/Documents/github/muller_diagrams/tests/data/tables_ggmuller/m5_correct.ggmuller.populations.tsv"
filename <- args[1]
start_positions <- 0.5
pop_df <- read.csv(filename, sep = "\t", header = TRUE)
pop_df <- add_start_points(pop_df)

write.table(pop_df, args[2], sep = "\t", row.names = FALSE)
