
library(ggplot2)
library(ggmuller)

population <- read.table("/home/cld100/Documents/github/muller_diagrams/example/B1_muller_try1.ggmuller.populations.tsv", header=TRUE)
edges <- read.table("/home/cld100/Documents/github/muller_diagrams/example/B1_muller_try1.ggmuller.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
write.csv(Muller_df, "example/supplementary_files/B1_muller_try1.muller.tsv", sep = "\t", col.names = NA)
palette <- c("#333333","#e6194b","#3cb44b","#ffe119","#4363d8","#f58231","#911eb4","#46f0f0","#f032e6","#bcf60c")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
geom_area() +
theme(legend.position = "right") +
guides(linetype = FALSE, color = FALSE) +
scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
scale_x_continuous(expand=c(0,0)) +
scale_fill_manual(name = "Identity", values = palette) +
scale_color_manual(values = palette)


ggsave("/home/cld100/Documents/github/muller_diagrams/example/B1_muller_try1.muller.basic.png", height = 10, width = 10)
