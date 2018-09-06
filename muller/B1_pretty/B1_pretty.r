
library("ggplot2")
library("ggmuller")

population <- read.table("/home/cld100/Documents/github/muller_diagrams/muller/B1_pretty/B1_pretty.ggmuller_populations.csv", header=TRUE)
edges <- read.table("/home/cld100/Documents/github/muller_diagrams/muller/B1_pretty/B1_pretty.ggmuller_edges.csv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
Muller_plot(Muller_df)

ggsave("/home/cld100/Documents/github/muller_diagrams/muller/B1_pretty/B1_pretty.muller.png", height = 10, width = 10)

