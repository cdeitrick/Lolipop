library("dplyr")
find_start_node <- function(edges) {
  start <- sort(edges$Parent)[1] # reasonable guess
  if(is.factor(start)) start <- levels(start)[start]
  repeat {
    if(move_up(edges, start) == start) break
    start <- move_up(edges, start)
  }
  return(start)
}
move_up <- function(edges, identity) {
  if(!(identity %in% edges$Identity) & !(identity %in% edges$Parent)) stop("Invalid identity.")
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

path_vector <- function(edges) {
  n <- 1
  path <- find_start_node(edges)
  debugvaluepath = path[n]
  upped <- FALSE
  repeat {
    if(!upped) repeat { # a downwards move should never follow an upwards move
      n <- n + 1
      path[n] <- move_down(edges, path[n - 1])
      upped <- FALSE
      if(path[n] == path[n - 1]) break
    }
    if(move_right(edges, path[n]) != path[n]) {
      n <- n + 1
      path[n] <- move_right(edges, path[n - 1])
      upped <- FALSE
    } else if(move_up(edges, path[n]) != path[n]) {
      n <- n + 1
      path[n] <- move_up(edges, path[n - 1])
      upped <- TRUE
    }
    if(path[n] == path[1]) break
    if(n > 2 * dim(edges)[1] + 2) stop("Error: stuck in a loop")
    if(max(table(path) > 2)) stop("Error: adjacency matrix seems to include loops.")
  }
  if(length(path) != 2 * dim(edges)[1] + 2) stop("Error: adjacency matrix seems to be bipartite.")
  return(path)
}

args = commandArgs(trailingOnly=TRUE)

table_edges <- read.table(args[1], sep = "\t", header = TRUE)

result <- path_vector(table_edges)
write(result, file = args[2])

