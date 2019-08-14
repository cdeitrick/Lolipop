library("dplyr")


find_start_positions <- function(pop_df, start_positions = 0.5) {
	# set small time interval:
	all_gens_list <- unique(pop_df$Generation)
	delta <- abs(min(1E-2 * min(diff(all_gens_list)), 1E-4 * (max(all_gens_list) - min(all_gens_list))))
	start_positions <- max(start_positions, delta)
	start_positions <- min(start_positions, 1 - delta)
	
	return(start_positions)
}

args = commandArgs(trailingOnly=TRUE)
#filename <- "/home/cld100/Documents/github/muller_diagrams/5G/tables/5_genotypes.timeseries.ggmuller.populations.tsv"
filename <- args[1]
pop_df <- read.csv(filename, sep = "\t", header = TRUE)
result <- find_start_positions(pop_df)

write(result, args[2])