working_dir <- "C:/Users/a499a400/Dropbox/chan"
charset <- "Char.txt"
tree <- "ExaBayes_ConsensusExtendedMajorityRuleNexus.contree.tre"

#SNPs_to_clades <- function(working_dir,charset,tree) {

library(stringr)
setwd(working_dir)
charsetable <- as.matrix(read.table(charset,sep="\t"))
treeset <- readLines(tree)

x <- 0

temptaxa <- NULL
temptree <- NULL

for (i in 1:(length(treeset))) {
if (grepl("taxlabels",treeset[i],fixed=TRUE)==TRUE) {
x <- 1
}
if (grepl(";",treeset[i],fixed=TRUE)==TRUE) {
x <- 0
}
if (grepl("(",treeset[i],fixed=TRUE)==TRUE) {
x <- 2
}
if (x==1) {
temptaxa <- rbind(temptaxa,treeset[i])
}
if (x==2) {
temptree <- rbind(temptree,treeset[i])
}
}

temptaxa <- temptaxa[-1,1]
subsample <- seq(3,(length(temptaxa)*3),3)
temptaxa <- unlist(strsplit(temptaxa,"\t",fixed=TRUE))[subsample]

temptree <- unlist(strsplit(temptree," ",fixed=TRUE))
temptree <- temptree[length(temptree)]
temptree <- gsub("\\[.*?\\]","",temptree,fixed=FALSE)
temptree <- gsub("\\:[0-9]+.[0-9]+[e-]*[0-9]*","",temptree,fixed=FALSE)
temptree <- gsub(";","",temptree,fixed=TRUE)
temptree <- unlist(strsplit(temptree,","))

for (i in 1:(length(temptaxa))) {
for (j in 1:(length(temptree))) {
temppattern <- gsub("\\(","",temptree[j])
temppattern <- gsub("\\)","",temppattern)
if(temppattern==i) {
temptree[j] <- gsub(i,temptaxa[i],temptree[j])
break
}
}
}

nodenames <- c(temptaxa,charsetable[1,-1])
nodenames <- nodenames[duplicated(nodenames)==FALSE]
nodenames <- nodenames[(length(temptaxa)+1):(length(nodenames))]

nodenames <- sort(as.numeric(nodenames))

recordtaxa <- matrix(NA, ncol=2,nrow=(length(nodenames)))
recordtaxa[,1] <- nodenames

m <- 0

for (i in 1:(length(temptree))) {
numopens <- nchar(temptree[i])-nchar(gsub("(","",temptree[i],fixed=TRUE))
if (numopens > 0 ) {
for (k in 1:numopens) {
x <- 1+numopens-k
m <- m + 1
j <- i + 1
while (j <= (length(temptree))) {
x <- x+nchar(temptree[j])-nchar(gsub("(","",temptree[j],fixed=TRUE))-nchar(temptree[j])+nchar(gsub(")","",temptree[j],fixed=TRUE))
if (x<=0) {
break
}
j <- j+1
}
torecord <- paste(temptree[i:j],collapse=",")
recordtaxa[m,2] <- torecord
}
}
}


recordtaxa[1:dim(recordtaxa)[1],2] <- gsub("(","",recordtaxa[1:dim(recordtaxa)[1],2],fixed=TRUE)
recordtaxa[1:dim(recordtaxa)[1],2] <- gsub(")","",recordtaxa[1:dim(recordtaxa)[1],2],fixed=TRUE)

recordtaxa <- recordtaxa[order(nchar(recordtaxa[,2])),]

for (i in 1:(dim(recordtaxa)[1])) {
for (j in (i+1):(dim(recordtaxa)[1])) {
recordtaxa[j,2] <- gsub(recordtaxa[i,2],recordtaxa[i,1],recordtaxa[j,2],fixed=TRUE)
}
}

write.table(recordtaxa,"list_of_clades_in_tree.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)


