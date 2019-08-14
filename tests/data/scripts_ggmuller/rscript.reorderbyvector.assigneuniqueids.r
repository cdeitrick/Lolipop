
# add unique id column to the vector:

generateUniqueIdsGenotype <- function(vector) {
	dup <- duplicated(vector)
	vector <- as.data.frame(vector)
	vector$count <- 1:nrow(vector)
	B <- data.frame(Unique_id = apply(vector, 1, function(x)
		if(dup[as.numeric(x["count"])]) paste0(x["vector"], "a") 
		else x["vector"]))
	B <- data.frame(lapply(B, as.character), stringsAsFactors=FALSE)
	vector <- cbind(vector[-ncol(vector)], B)
}

args = commandArgs(trailingOnly=TRUE)

vectorstrings <- args[1]
#vectorstrings <- "genotype-0,genotype-1,genotype2"
vector <- strsplit(vectorstrings, ",")[[1]]

result <- generateUniqueIdsGenotype(vector)

write.table(result, file = args[2], sep = "\t")
