library("dplyr")

cleanEdges <- function(edges, lookup) {
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
}

args = commandArgs(trailingOnly=TRUE)

edges_table <- read.table(args[1], sep = "\t", header = TRUE)
lookup_table <- read.table(args[2], sep = "\t", header = TRUE)

clean_edges <- cleanEdges(edges_table, lookup_table)

write.table(clean_edges, args[3], sep = "\t")
