
library("ggplot2")
library("ggmuller")

population <- read.table("/home/cld100/Documents/github/muller_diagrams/muller/example/example.ggmuller_populations.tsv", header=TRUE)
edges <- read.table("/home/cld100/Documents/github/muller_diagrams/muller/example/example.ggmuller_edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
Muller_plot(Muller_df)

ggsave("/home/cld100/Documents/github/muller_diagrams/muller/example/example.muller.png", height = 10, width = 10)

