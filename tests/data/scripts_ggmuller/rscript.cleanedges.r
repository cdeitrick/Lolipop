

cleanEdges <- function(pop_df, lookup) {
	# replace each genotype name in adjacency matrix with corresponding Age:
	edges <- filter_(edges, ~Identity %in% lookup$Identity)
	edges <- left_join(edges, lookup, by = "Identity")
	edges <- select_(edges, ~-Identity)
	colnames(edges) <- c("Parent", "Identity")
	edges <- arrange_(edges, ~Identity)
	colnames(lookup)[1] <- "Parent"
	edges <- left_join(edges, lookup, by = "Parent")
	edges$Parent <- edges$Age
	edges <- select_(edges, ~-Age)
	return(edges)
}

args = commandArgs(trailingOnly=TRUE)

table_edges <- read.table(args[1], sep = "\t", header = "TRUE")
table_lookup <- read.table(args[2], sep = "\t", header = TRUE)

cleaned_table <- cleanEdges(table_edges, table_lookup)

write.table(cleaned_table, args[3], sep = "\t")