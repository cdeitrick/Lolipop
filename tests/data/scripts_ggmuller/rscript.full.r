cleanPopulationData <- function(pop_df) {
	# add missing population values:
	if(dim(pop_df)[1] != length(unique(pop_df$Identity)) * length(unique(pop_df$Generation))) {
		added_rows <- expand.grid(Identity = unique(pop_df$Identity), Generation = unique(pop_df$Generation))
		added_props <- group_by(pop_df, Identity) %>% 
			slice(1) %>% 
			ungroup() %>% 
			select(-one_of("Generation", "Population"))
		added_rows <- merge(added_rows, added_props, all = TRUE)
		pop_df <- merge(added_rows, pop_df, all = TRUE)
		pop_df[is.na(pop_df$Population), "Population"] <- 0
		pop_df <- arrange_(pop_df, ~Generation)
		warning("missing population sizes replaced by zeroes")
	}
	
	if (!missing(add_zeroes)) {
		warning("argument add_zeroes is deprecated (it is now always TRUE).", 
				call. = FALSE)
	}
	if (!missing(smooth_start_points)) {
		warning("argument smooth_start_points is deprecated (it is now always TRUE).", 
				call. = FALSE)
	}
	if (!missing(threshold)) {
		warning("argument threshold is deprecated (use cutoff instead, noting that genotypes whose abundance never exceeds the cutoff value are removed, 
            whereas previously genotypes whose abundance never exceeded *twice* the threshold value were removed).", 
				call. = FALSE)
		if (missing(cutoff)) cutoff <- threshold * 2
	}
	
	# check/set column names:
	if(!("Generation" %in% colnames(pop_df)) | !("Identity" %in% colnames(pop_df)) | !("Generation" %in% colnames(pop_df))) 
		stop("colnames(pop_df) must contain Generation (or Time), Identity and Population")
	
	if(!is.na(edges)[1]) {
		set1 <- unique(pop_df$Identity)
		set2 <- unique(edges$Identity)
		set3 <- unique(edges$Parent)
		# check that pop_df and edges have compatible Identity values:
		if(length(setdiff(set1, set2)) != 1) stop("Identity values in edges must match Identity values in pop_df, excluding the original genotype (which has no parent)")
		# check that Parent and Identity values in edges are consistent:
		if(length(setdiff(set3, set2)) != 1) stop("Parent values in edges must also appear as Identity values in edges, excluding the original genotype (which has no parent)")
	}
	
	if(!is.na(edges)[1]) {
		if("phylo" %in% class(edges)) {
			collapse.singles(edges)
			edges <- edges$edge
		}
		edges <- na.omit(edges) # remove any rows containing NA
		colnames(edges) <- c("Parent", "Identity")
		if(is.factor(edges$Parent)) edges$Parent <- levels(edges$Parent)[edges$Parent]
		if(is.factor(edges$Identity)) edges$Identity <- levels(edges$Identity)[edges$Identity]
	}
	
	return(pop_df)
}

find_start_positions <- function(pop_df, start_positions = 0.5) {
	# set small time interval:
	all_gens_list <- unique(pop_df$Generation)
	delta <- abs(min(1E-2 * min(diff(all_gens_list)), 1E-4 * (max(all_gens_list) - min(all_gens_list))))
	start_positions <- max(start_positions, delta)
	start_positions <- min(start_positions, 1 - delta)
	
	return(start_positions)
}

filename <- "/home/cld100/Documents/github/muller_diagrams/5G/tables/5_genotypes.timeseries.ggmuller.populations.tsv"
pop_df <- read.csv(filename, sep = "\t", header = TRUE)
pop_df <- cleanPopulationData(pop_df)
result <- find_start_positions(pop_df)
