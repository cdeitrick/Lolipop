library("dplyr")

move_up <- function(edges, identity) {
  if(!(identity %in% edges$Identity) & !(identity %in% edges$Parent)) stop("Invalid identity.")
  parent <- filter_(edges, ~Identity == identity)$Parent
  if(length(parent) == 0) return(identity) # if it is the initial genotype then don't move
  if(is.factor(parent)) parent <- levels(parent)[parent]
  return(parent)
}

find_start_node <- function(edges) {
  start <- sort(edges$Parent)[1] # reasonable guess
  if(is.factor(start)) start <- levels(start)[start]
  repeat {
    if(move_up(edges, start) == start) break
    start <- move_up(edges, start)
  }
  return(start)
}

args = commandArgs(trailingOnly=TRUE)

table_edges <- read.table(args[1], sep = "\t", header = TRUE)
print(table_edges)
result <- find_start_node(table_edges)

write(result, file = args[2])
