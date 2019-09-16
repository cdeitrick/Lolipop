
library(ggplot2)
library(ggmuller)

population <- read.table("/home/cld100/Documents/github/muller_diagrams/strongselection/tables/model.strongselection.genotype.ggmuller.populations.tsv", header=TRUE)
edges <- read.table("/home/cld100/Documents/github/muller_diagrams/strongselection/tables/model.strongselection.genotype.ggmuller.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
palette <- c("#FFFFFF","#c7e9c0","#6aaed6","#e32f27","#228a44","#fca082","#73c476")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
geom_area() +
theme(legend.position = "right") +
guides(linetype = FALSE, color = FALSE) +
scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
scale_x_continuous(expand=c(0,0)) +
scale_fill_manual(name = "Identity", values = palette) +
scale_color_manual(values = palette)

ggsave("/home/cld100/Documents/github/muller_diagrams/strongselection/graphics/model.strongselection.genotype.rscript.png", height = 10, width = 10)
