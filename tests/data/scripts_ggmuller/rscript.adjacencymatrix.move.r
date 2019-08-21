library("dplyr")

move_up <- function(edges, identity) {
	if(!(identity %in% edges$Identity) & !(identity %in% edges$Parent)) stop(paste("Invalid identity: ", identity))
	parent <- filter_(edges, ~Identity == identity)$Parent
	if(length(parent) == 0) return(identity) # if it is the initial genotype then don't move
	if(is.factor(parent)) parent <- levels(parent)[parent]
	return(parent)
}

move_right <- function(edges, identity) {
	if(!(identity %in% edges$Identity) & !(identity %in% edges$Parent)) stop("Invalid identity.")
	parent <- filter_(edges, ~Identity == identity)$Parent
	if(length(parent) == 0) return(identity) # if it is the initial genotype then don't move
	siblings <- sort(filter_(edges, ~Parent == parent)$Identity)
	siblings <- siblings[which(siblings == identity) + 1]
	if(length(siblings) == 0) return(identity) # if it is the initial genotype then don't move
	if(is.na(siblings)) return(identity) # if it is the initial genotype then don't move
	if(is.factor(siblings)) siblings <- levels(siblings)[siblings]
	return(siblings)
}


move_down <- function(edges, parent) {
	if(!(parent %in% edges$Identity) & !(parent %in% edges$Parent)) stop("Invalid parent.")
	daughters <- filter_(edges, ~Parent == parent)$Identity
	if(length(daughters) == 0) return(parent) # if it is not a parent then don't move
	if(is.factor(daughters)) daughters <- levels(daughters)[daughters]
	return(sort(daughters)[1])
}

debugvalue <- FALSE

args = commandArgs(trailingOnly=TRUE)


filename <- ifelse(debugvalue, "/tmp/pytest-of-cld100/pytest-170/test_moves_B1_muller_try1_ggmu0/adjacencymatrix", args[1]) 
edges <- read.table(filename, sep = "\t", header = TRUE)
test_value <- ifelse(debugvalue, 2, as.numeric(args[2]))

resultup <- move_up(edges, test_value)
resultdown <- move_down(edges, test_value)
resultright <- move_right(edges, test_value)

result <- c(resultup, resultdown, resultright)

if (!debugvalue) {
	write(result, args[3])
}
#

