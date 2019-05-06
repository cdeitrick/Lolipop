
library(ggplot2)
library(ggmuller)

population <- read.table("/home/cld100/Documents/github/muller_diagrams/R3/tables/example.ggmuller.populations.tsv", header=TRUE)
edges <- read.table("/home/cld100/Documents/github/muller_diagrams/R3/tables/example.ggmuller.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
write.csv(Muller_df, "/home/cld100/Documents/github/muller_diagrams/R3/tables/example.muller.csv", sep = "\t", col.names = NA)
palette <- c("#FFFFFF","#6aaed6","#540000","#ffffc8","#fcbba1","#ffff5a","#2f0000","#5b3495","#ff3400","#ffa500","#c6dbef","#aedea7","#2070b4","#9e0000","#ffff91","#ff5c00","#ffca00","#796eb2","#c6c7e1","#ca181d","#ff8000","#c50000","#37a055","#ff1000","#9e9ac8","#fb694a","#790000","#ea0000","#e8e7f2","#ffff22","#ffef00","#959595")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
geom_area() +
theme(legend.position = "right") +
guides(linetype = FALSE, color = FALSE) +
scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
scale_x_continuous(expand=c(0,0)) +
scale_fill_manual(name = "Identity", values = palette) +
scale_color_manual(values = palette)

ggsave("/home/cld100/Documents/github/muller_diagrams/R3/graphics/clade/example.muller.basic.png", height = 10, width = 10)
