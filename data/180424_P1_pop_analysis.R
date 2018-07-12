#this is the P1 population of the LTE that I am trying to automate the analysis of

#The input for this script is a csv file of a breseq output. Can use the breseq_cat script that nate/Chris Dietrick made to get an excel file. You do have to do a little bit of manipulation on this file in excel first, described in comments below. 

#this script will take the list of mutations and subtract mutations:1. found in the ancestor, 2. that don't reach a cumulative frequency of 10%, 3. that only appear at one time point, 4. that are fixed from the first measured time point, 5. that are almost fixed from the first time point, and 6. that do not change in frequency by at least 10% over the course of the experiment.

#the output of this script is an excel sheet that is formatted to go directly into Katya's matlab scripts to define mutational cohorts, infer ancestry, and make muller plots. there are many places that data frames can be printed for other purposes. 

library("vegan")
library("plyr")
library("RColorBrewer")
library("ggplot2")
library("data.table")
library("dplyr")
library("reshape2")
library("xlsx")
library("scales")

theme_set(theme_bw())

#first thing is downloading time course breseq data from beagle
#/home/kah231/scripts/BreseqCat.py -d /home/kah231/PA/analysis/LTEpitt_seq/P1
#and the ancester
#/home/kah231/scripts/SingleBreseqCat.py -f /home/kah231/PA/analysis/LTEpitt_seq/KBH5_WT/output/index.html
#
#
##get rid of spaces in all columns except description.
#convert arrows ← and → o < and > get rid of all commas, remove Â, ‑,–
#need to make sure days are just days not "P117" etc.
##saved the snp tab as a csv format


#FIle locations
#/Users/katrina/Desktop/working
#called:
#  P1_Breseq_Output.csv
#  KBH5_WT_Breseq_Output.csv

setwd("/Users/katrina/Desktop/working")

#ancestral/background SNPs
ancestor_snps <- read.csv("KBH5_WT_Breseq_Output.csv", header=TRUE)
head(ancestor_snps)
View(ancestor_snps) #want SeqID column because that is where the positions are

#P1 population SNPs
P1_snps <- read.csv("P1_Breseq_Output.csv",header=TRUE)
View(P1_snps) #again SeqID is what we want to match
nrow(P1_snps) #4494 total

#get the SNPs that are found in both files according to the position numbers that are
#actually stored in the SeqID column
#filter1 = no mutations found in the ancestral clone
P1_filter1 <- P1_snps[ !(P1_snps$SeqID %in% ancestor_snps$SeqID), ]
#see how many rows are in the data frame now
nrow(P1_filter1) #3534
#create a data frame of the shared mutations (the ones taken out with filter 1)
P1_filter1_refmatch <- P1_snps[ (P1_snps$SeqID %in% ancestor_snps$SeqID), ]
nrow(P1_filter1_refmatch)#960
write.csv(P1_filter1_refmatch, file = "P1_ancestral.csv")


#will create a df with taking out based on the gene name. I don't think this will be what I want, I think it will be too strict but I am doing it anyway
#I don't like using this because some genes really can get more than one mutation, especially if it is a big gene
P1_genelevelfilter <- P1_snps[ !(P1_snps$Gene %in% ancestor_snps$Gene), ]
nrow(P1_genelevelfilter)#3272
#now see how many it took out
P1_genelevelfilter_refmatch <- P1_snps[ (P1_snps$Gene %in% ancestor_snps$Gene), ]
nrow(P1_genelevelfilter_refmatch)#1222

#this is done to the original data set through filter1
#remove '%' symbol
P1_filter1$Mutation <- gsub( "%", "", as.character(P1_filter1$Mutation), n)
P1_filter1$Mutation <- as.numeric(as.character(P1_filter1$Mutation))

#combine annotation::gene::description in one column
P1_filter1$desc_gene_annot <-  paste(P1_filter1$Annotation,P1_filter1$Gene,P1_filter1$Description, sep="::")
P1_filter1$desc_gene_annot <- as.factor(P1_filter1$desc_gene_annot)
#combine desc_gene_annot and position (SeqID) in one column
P1_filter1$details <- paste(P1_filter1$desc_gene_annot, P1_filter1$SeqID, sep=";;")
P1_filter1$details <- as.factor(P1_filter1$details)
#combine details with the actual mutation
P1_filter1$info <- paste(P1_filter1$details, P1_filter1$Position, sep=";;")
P1_filter1$info <- as.factor(P1_filter1$info)
View(P1_filter1)

#melt data frame for casting
m_P1_filter1 <- melt(P1_filter1, id=c("Sample","Evidence","SeqID","Position","Annotation","Gene","Description","desc_gene_annot", "details", "info"),measure.vars = c("Mutation"))
head(m_P1_filter1)
View(m_P1_filter1)

#cast data frame - organizing with each mutation as the rows and the frequency of that mutation on a given day as the columns
P1_cast1 <- t(dcast(m_P1_filter1,Sample~info,mean, value.var = "value",fill=0))
P1_cast1 <- as.data.frame(P1_cast1,header=TRUE)
colnames(P1_cast1) <- as.character(unlist(P1_cast1[1,]))
colnames(P1_cast1) 
View(P1_cast1)
P1_cast1$"0" <- 0.0
View(P1_cast1)
P1_cast1 <- P1_cast1[-1,]
View(P1_cast1)

#need to reorder the columns in ascending order
P1_column_order <- c("0","17", "44","66","90")
colnames(P1_cast1)
setcolorder(P1_cast1, P1_column_order)
View(P1_cast1)

nrow(P1_cast1)#1946

#transpose the matrix
t_P1 <- as.data.frame(t(P1_cast1))
View(t_P1)
#figure out what class the frequency values are in the matrix - they need to be numeric
class(t_P1[2,2]) #factor
ncol(t_P1)#1946 <- for sanity this check should match up with what you found previously for a nrow count

#convert frequency values to numeric class - start as "character"
t_P1[,2:1946] <-(apply(t_P1[,2:1946], 2, function(x) as.numeric(as.character(x)))) 
P1 <- transpose(t_P1[,2:1946])
colnames(P1) <- rownames(t_P1)
rownames(P1) <- colnames(t_P1[,2:1946])
View(P1)

class(P1[2,2]) #yay, it's numeric now!!

#adds a count number on each row (mutation) that tells how many columns (days) have a frequency above 0
P1$count <- rowSums(P1!=0.0)
#sums up the rows to tell you the total % that you get - flaw is that it also adds the count, so need to subtract that value
P1$Sums <- rowSums(P1)-P1[,6]
nrow(P1)#1945
View(P1) #these last two are mostly just sanity checks

#filter 2!!!!! select only rows with greater than 10% total frequency
P1_filter2 <- (subset(P1, P1$Sums >= 10)) #greater than 10%
nrow(P1_filter2)#1028

P1_filter2_out <- (subset(P1, P1$Sums < 10)) #put all of the filtered out mutations in one place

#filter 3!! select only rows that appear in more than 1 day
P1_filter3 <- (subset(P1_filter2, P1_filter2$count > 1))
nrow(P1_filter3)#826
P1_filter3_out <- (subset(P1_filter2, P1_filter2$count <= 1))
View(P1_filter3)

#filter 4 -- remove all mutations at 100% across all measured time points
#4 time points so 400 value --> problem with this is mutations can start at 100 and dip slightly
P1_filter4 <- (subset(P1_filter3, P1_filter3$Sums < 400))
nrow(P1_filter4) #816
P1_filter4_out <- (subset(P1_filter3, P1_filter3$Sums >= 400))

#filter out if the first time points that are 95% or above
P1_filter5 <- (subset(P1_filter4, P1_filter4$"17" < 95))
nrow(P1_filter5) #814
P1_filter5_out <- (subset(P1_filter4, P1_filter4$"17" >= 95))

#filter out if the HIGHEST frequency isn't 10, not if the combined total frequency doesn't get to 10

#filter out if the change in frequency from the first time point to the last time point does not change by at least 10%, having the additive value be above 10 isn't stringent enough.
############need to change this because it does not do what I thought it did
#P1_filter6 <- (subset(P1_filter5, (P1_filter5$"17"+P1_filter5$"90")/2 >= 10))
#nrow(P1_filter6)#248
#filter6_removed <- (subset(P1_filter5, (P1_filter5$"17"+P1_filter5$"90")/2 < 10))
#View(filter6_removed)
#I really should be checking what mutations are being taken out with each filter... 
#View(P1_filter6)
###############
#what about subtracting 90 from 17 and thenn if the absolute value isn't above 10 filter it out
P1_filter62 <- (subset(P1_filter5, abs(P1_filter5$"17"-P1_filter5$"90") >= 10))
nrow(P1_filter62) #104
P1_filter62_out <- (subset(P1_filter5, abs(P1_filter5$"17"-P1_filter5$"90") < 10))
View(P1_filter62)

not_real_mutations <- rbind(P1_filter2_out, P1_filter3_out, P1_filter4_out, P1_filter5_out, P1_filter62_out)
write.csv(not_real_mutations, file = "P1_Filtered_out.csv")


ncol(P1_filter62)#7
P1_filter7 <- P1_filter62[,-c(6,7)] #remove columns with count and sums
View(P1_filter7)

#split the row names into columns again
P1_split <- P1_filter7
P1_split$info <- rownames(P1_filter7)
View(P1_split)

P1_split4 = transform(P1_split, info =colsplit(P1_split$info,';;', names = c('desc_gene_annot','position', 'Mutation')))
#View(P1_split4)
P1_split4 <- as.data.frame(P1_split4)
colnames(P1_split4)#6
info <- (P1_split4$info) 
head(info)
P1_split4 <- P1_split4[,c(1:5)]
View(P1_split4)

#P1_split4$desc_gene_annot <- rownames(P1_split4)
#P1_split4$desc_gene_annot <- gsub(";;.*","",P1_split4$desc_gene_annot)
#does_this_work <- merge(P1_split4, info, by="desc_gene_annot")
#this does work but I like my way better because I know what it is doing and I still have the rownames
#View(does_this_work )

View(P1_split4)
View(info)
P1_split5 <- P1_split4
P1_split5$Position <- info$position
P1_split5$Mutation <- info$Mutation
View(P1_split5)
#P1_split6 <- P1_split5[,-6]
#View(P1_split6)

#rename columns after splitting
#colnames(P1_split6)
colnames(P1_split5) <- c("0","17", "44","66","90","Position", "Mutation")
View(P1_split5)

#write this to a file that I can find
write.csv(P1_split5,file="P1_allfilters.csv")

#oh plotting....
#have to melt it first...
#melt - with whatever i want to keep 
#transform
View(P1_split5)
nrow(P1_split5)#92

t_P1_plot <- t(P1_split5)
nrow(t_P1_plot)#7
t_P1_plot_2 <- t_P1_plot[-c(6,7),]
m_P1_plot <- melt(t_P1_plot_2, id=c("desc_gene_annot"),Value.name = "frequency")
head(m_P1_plot)

#plot
#plot(NA, xlim=c(0,100), ylim=c(0,100))
#lines(m_P1_plot$X1, m_P1_plot$value)
#lines(c(0,17,44,66,90), t_final_P1_filter3[,1], add=TRUE)
#lines(c(0,17,44,66,90), t_final_P1_filter3[,2], add=TRUE)
#this is how I could add things one at a time

#but I want to know how to use ggplot

colnames(m_P1_plot) <- c("day", "mutation","value")
m_P1_plot$value <- as.numeric(as.character(m_P1_plot$value)) # you have to change from a factor to a number this way. You have to co to a character first always. if I had just gone to a number it would have given me the level of factor that the previous data point was. 
View(m_P1_plot)

ggplot(m_P1_plot,aes(x=day,y=value,color=mutation)) +theme(text = element_text(size=20),legend.text=element_text(size=10),legend.position="none") +geom_point(size=4) +geom_line()


#this data set is going to be formatted to go into katya's matlab scripts. The data frame needs specific columns with particular data. they will have to be in a specific order and named correctly, but first I just need to make them holding the correct information. 
View(P1_split5)
P1_Muller <- P1_split5[,c(1:5)]
P1_Muller$Population <- "P1"
P1_Muller$Population2 <- 1 #make sure this is a number
P1_Muller$Chromosome <- 1 #make sure this is a number
P1_Muller$Position <- P1_split5$Position
P1_Muller$Class <- "SNP"
P1_Muller$Mutation <- P1_split5$Mutation
P1_Muller$Gene <- ""
P1_Muller$AminoAcid <- ""
P1_Muller$Class2 <- ""
P1_Muller$Amino <- ""
P1_Muller$NearestDownstreamGene <- ""
P1_Muller$Distance <- ""
P1_Muller$Trajectory <- 1:nrow(P1_Muller)

colnames(P1_Muller)
#now put the columns in the correct order
Muller_col_order <- c("Population", "Population2","Trajectory","Chromosome","Position","Class","Mutation","Gene","AminoAcid","Class2","Amino","NearestDownstreamGene","Distance","0","17","44","66","90")
setcolorder(P1_Muller,Muller_col_order)

#now I need to name them what they are actually supposed to be named
colnames(P1_Muller)

colnames(P1_Muller) <-c("Population", "Population number","Trajectory","Chromosome","Position","Class","Mutation","Gene","Amino Acid","Class","Amino","Nearest Downstream Gene","Distance","0","17","44","66","90")

View(P1_Muller)
#need to remove the rownames
#rownames(P1_Muller) <- c()
#View(P1_Muller)
#decided not to do this because I can just print without including the row names. this will also just print out the row names, they will just be the numbers instead of the descriptions


#latest problem is that the frequencies need to be percentages and the column names for the frequencies need to be numbers. 
#first solve the frequencies to percentages problem - should be able to do with scales package 
#would like to keep 2 decimal points if possible

P1_Muller_try <- P1_Muller
P1_Muller_try$`90` <- as.numeric(as.character(P1_Muller_try$`90`))
P1_Muller_try$`66` <- as.numeric(as.character(P1_Muller_try$`66`))
P1_Muller_try$`44` <- as.numeric(as.character(P1_Muller_try$`44`))
P1_Muller_try$`17` <- as.numeric(as.character(P1_Muller_try$`17`))
P1_Muller_try$`0` <- as.numeric(as.character(P1_Muller_try$`0`))

#the following does work. I divide everything by 100 so that when in excel I can change to "percent" type and it will be the correct value. Remember to keep 1 decimal when changing the type in excel or it will round everything. 
P1_Muller_try$`90` <- (P1_Muller_try$`90`/100)
P1_Muller_try$`66` <- (P1_Muller_try$`66`/100)
P1_Muller_try$`44` <- (P1_Muller_try$`44`/100)
P1_Muller_try$`17` <- (P1_Muller_try$`17`/100)
P1_Muller_try$`0` <- (P1_Muller_try$`0`/100)


View(P1_Muller_try)

#now to write this file so that I can use it as an input to matlab. The matlab file requires it to be a .xlsx file so I can just write to that type of file. need to make sure that I don't print out the row names or they will be the first column. I do NEED the column names though. 
write.xlsx(P1_Muller_try, file="P1_Muller.xlsx", row.names = FALSE)
#this file needs to be loaded into matlab for Katya's scripts.






#see what the last filter took out to make sure the mutations are what I want them to be. This check needs to be done on all mutations that are taken out from step 1 on. 
filter62_removed <- (subset(P1_filter5, abs(P1_filter5$"17"-P1_filter5$"90") < 10))
View(filter62_removed)#yay it removed what I wanted it to!!
nrow(filter62_removed)#587
#I need to look at the pileups for these mutations to se why they are coming up in my populations. 
#there are way more of these "nonsense" mutations than there dhould be if it were just normal sequencing errors. so the question is where is the source of the variation. I will look at the pileup at the areas of these mutations in IGV to try to determine the answer to this. 

split_removed <- filter62_removed
split_removed$info <- rownames(filter62_removed)
split_removed2 = transform(split_removed, info =colsplit(split_removed$info,';;', names = c('desc_gene_annot','position', 'Mutation')))
split_removed2 <- as.data.frame(split_removed2)
colnames(split_removed2)
View(split_removed2)
info2 <-(split_removed2$info)
View(info2)
ncol(split_removed2)
split_removed2 <- split_removed2[,-c(8)]
View(split_removed2)
colnames(info2)
split_removed2$position <- info2$position
split_removed2$Mutation <- info2$Mutation
View(split_removed2)
write.csv(split_removed2, file = "P1_removed.csv")


#######################################
#the following is a section of MatLab code (so won't work to run in R)

#problem I ran into is that there is apparently something wrong with the excel file... it won't read in
#manipulations done to excel file once exported 
#change column names of frequencies to "number" type
#change the frequencies to "percentage" type
#
names = ["P1_muller.xlsx"]; %take in the excel file you got from R
sheets = ["Sheet1"]; %it is an excel workbook so you have to tell it which sheet it's on
tNum = 5; %the number of timepoints
timepoints = [0,17,44,66,90]; %the names of the timepoints
xUnits = ["Days"]; %units of the timepoints
time_series_import_KMF %this reads in the data
get_genotypes_KMF_KBH %this will determine the genotypes --> need to learn how it is defining a genotype
genotype_plots_KMF %this plots the genotypes. It will give 3 plots. the first is unclustered, there is an intermediate, and the third is the final fully clustered data set.


order_clusters_KMF % this will determine the order in which the genotypes showed up
ordered_cluster_plots_KMF %visualize the final clusters with the ancestry incorporated --> should save this figure.

frequencies = squeeze(genneststotal(1, any(squeeze(genneststotal(1, :, :)), 2), 2:trajSize)); % define the frequencies variable
nests = squeeze(genneststotal(1, any(squeeze(genneststotal(1, :, :)), 2), nameSize:end)); % define the nests variable
nests = nests(:, any(nests, 1)); %still defninng nests
muller_plots_KMF(frequencies, nests, timepoints) % problem with this is that it always outputs with the x axis saying "time in generations" I want to be able to change this to whatever I want to



