
library(ggplot2)
library(ggmuller)

population <- read.table("/home/cld100/Documents/github/example/tables/example.trajectories.ggmuller.populations.tsv", header=TRUE)
edges <- read.table("/home/cld100/Documents/github/example/tables/example.trajectories.ggmuller.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
write.csv(Muller_df, "/home/cld100/Documents/github/example/tables/example.trajectories.mullerdataframe.csv", sep = "\t", col.names = NA)
palette <- c("#FFFFFF","#4bb062","#3787c0","#2f974e","#ffca00","#fb694a","#b8e3b2","#98d594","#d3eecd","#73c476","#abd0e6","#006428","#157f3b","#9e9ac8","#ea0000","#e9f7e5")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
geom_area() +
theme(legend.position = "right") +
guides(linetype = FALSE, color = FALSE) +
scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
scale_x_continuous(expand=c(0,0)) +
scale_fill_manual(name = "Identity", values = palette) +
scale_color_manual(values = palette)

ggsave("/home/cld100/Documents/github/example/graphics/example.trajectories.rscript.png", height = 10, width = 10)
