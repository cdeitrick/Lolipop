
library(ggplot2)
library(ggmuller)

population <- read.table("/home/cld100/Documents/github/muller_diagrams/description/example/tables/example.ggmuller.populations.tsv", header=TRUE)
edges <- read.table("/home/cld100/Documents/github/muller_diagrams/description/example/tables/example.ggmuller.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
write.csv(Muller_df, "description/example/tables/example.muller.csv", sep = "\t", col.names = NA)
palette <- c("#FFFFFF","#334538","#87CA8C","#356A44","#4BAB63","#34573E","#E32F27","#368E4F","#357C49","#3787C1","#73C07F","#9AD49A","#37A155","#5EB570")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
geom_area() +
theme(legend.position = "right") +
guides(linetype = FALSE, color = FALSE) +
scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
scale_x_continuous(expand=c(0,0)) +
scale_fill_manual(name = "Identity", values = palette) +
scale_color_manual(values = palette)

ggsave("/home/cld100/Documents/github/muller_diagrams/description/example/graphics/clade/example.muller.basic.png", height = 10, width = 10)
