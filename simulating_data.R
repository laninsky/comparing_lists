library(stringr)
setwd("working_dir") # e.g. setwd("C:/Users/a499a400/Dropbox/Kaloula frogs")
input <- readLines("filename") # e.g. input <- readLines("froggies-25pctincomplete-filled-MAFFT-gblocks.nexus.phylip")
notaxa <- as.numeric(unlist(strsplit(input[1], " "))[1])

taxanames <- matrix(ncol=1,nrow=notaxa)

for (i in 1:notaxa) {
taxanames[i,1] <- unlist(strsplit(input[i+1]," "))[1]
}

sites <- 1000000
output <- matrix(ncol=sites,nrow=notaxa)
x <- c("C","T")

for (j in 1:sites) {
output[,j] <- sample(x,24,replace=TRUE)
}

tempmatrix <- matrix(ncol=2,nrow=notaxa)
tempmatrix[,1] <- taxanames[,1]

for (i in 1:notaxa) {
tempmatrix[i,2] <- paste(output[i,1:(dim(output)[2])],collapse="")
}

outputmatrix <- matrix(ncol=1,nrow=notaxa)
for (i in 1:notaxa) {
outputmatrix[i,1] <- paste(tempmatrix[i,],collapse=" ")
}

outputmatrix <- rbind("matrix",outputmatrix)
outputmatrix <- rbind("  format datatype=dna missing=? gap=-;",outputmatrix)
temp <- paste("  dimensions ntax=",notaxa," nchar=",sites,";",sep="")
outputmatrix <- rbind(temp,outputmatrix)
outputmatrix <- rbind("begin data;",outputmatrix)
outputmatrix <- rbind("#NEXUS",outputmatrix)
outputmatrix <- rbind(outputmatrix,";")
outputmatrix <- rbind(outputmatrix,"end;")

write.table(outputmatrix, "simulated.nexus",quote=FALSE, col.names=FALSE,row.names=FALSE)
