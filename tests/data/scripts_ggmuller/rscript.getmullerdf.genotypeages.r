library("dplyr")

get_genotype_ages <- function(pop_df) {
	# construct a dataframe with "Age" of each genotype:
	pop_df <- arrange_(pop_df, ~-Population)
	pop_df <- arrange_(pop_df, ~Generation)
	lookup <- group_by_(pop_df, ~Identity) %>% 
		filter_(~Population > 0 | Generation == max(Generation)) %>% 
		slice(1) %>% 
		arrange_(~Generation) %>% 
		ungroup()
	lookup <- mutate(lookup, Age = 1:dim(lookup)[1]) %>% 
		select_(~-c(Generation, Population))
	if(is.factor(lookup$Identity)) lookup$Identity <- levels(lookup$Identity)[lookup$Identity]
	lookup <- select_(lookup, ~c(Identity, Age))
	
	return(lookup)
}

args = commandArgs(trailingOnly=TRUE)
filename <- args[1]
pop_df <- read.table(filename, sep = "\t", header = TRUE)
					
# Add the sorting to match the python script
pop_df <- arrange_(pop_df, ~-Population)
pop_df <- arrange_(pop_df, ~Generation)

#write.table(pop_df, args[2], sep = "\t", row.names = FALSE)

#result <- add_start_points(pop_df)
lookup <- get_genotype_ages(pop_df)
#lookup <- pop_df
write.table(lookup, args[2], sep = "\t", row.names = FALSE)