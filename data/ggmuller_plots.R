library(ggmuller)
setwd("/Users/katrina/Desktop/working/ggmuller")

#edges <- read.csv("edges_try.csv", header=TRUE)
#View(edges)

#pop <- read.csv("pop_try.csv", header=TRUE)
#View(pop)

#Muller_try <- get_Muller_df(edges, pop, threshold = .005)
#Muller_plot(Muller_try)

#the above works as a test case so now I need to use my data from matLab to create these.

#my test is population P3, just because that was the data set that matlab had up when I started.
#to build the muller diagrams in matlab you need 3 data sets: frequencies, nests, and timepoints.
#gg muller combines these into two. one they call edges and one they call pops. edges denotes the ancestry of the lineages. this is matlab's nests, but simplified. Matlab stores the whole lineage whereas gg muller only requires the most recent ancestry. I believe that I am interpreting these correctly, but we will see when the figure pops out.

#p3_edges <- read.csv("edges_P3.csv", header = TRUE)
#View(p3_edges)

#p3_pop <- read.csv("pop_p3.csv", header= TRUE)
#View(p3_pop)
#p3_muller <- get_Muller_df(p3_edges, p3_pop, threshold = .0005)
#Muller_plot(p3_muller)


# p1 population
p1_edges <-  read.csv("edges_p1.csv", header = TRUE)
p1_pop <- read.csv("pop_p1.csv", header = TRUE)

View(p1_edges)
View(p1_pop)
p1_muller <- get_Muller_df(p1_edges, p1_pop, threshold= .005)
Muller_plot(p1_muller)

p1_pop2 <- read.csv("pop_p1_2.csv", header = TRUE)
View(p1_pop2)
p1_muller_2 <- get_Muller_df(p1_edges, p1_pop2, threshold = .005)
Muller_plot(p1_muller_2)

p1_pop3 <- read.csv("pop_p1_3.csv", header = TRUE)
View(p1_pop3)
p1_muller_3 <- get_Muller_df(p1_edges, p1_pop3, threshold = .005)
Muller_plot(p1_muller_3)

p1_edges2 <- read.csv("edges_p1_2.csv", header = TRUE)
p1_muller_4 <- get_Muller_df(p1_edges2, p1_pop2, threshold = .005)
Muller_plot(p1_muller_4)

p1_muller_5 <- get_Muller_df(p1_edges2, p1_pop, threshold = .005)
Muller_plot(p1_muller_5)
