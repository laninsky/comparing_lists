working_dir <- "C:/Users/a499a400/Dropbox/Kaloula frogs"
charset <- "origchar.txt"
tree <- "annotated_14July2015_concatenated.tre"

library(stringr)
setwd(working_dir)
charsetable <- as.matrix(read.table(charset,sep="\t"))
treeset <- readLines(tree)

#setting some variables for the loop below
x <- 0
temptaxa <- NULL
temptree <- NULL

#extracting the tree and taxa block from the treefile
for (i in 1:(length(treeset))) {
if (grepl("axlabels",treeset[i],fixed=TRUE)==TRUE) {
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

# making the tree and taxa variables prettier
temptaxa <- temptaxa[-1,1]
temptaxa <- gsub("\t"," ",temptaxa,fixed=TRUE)
temptaxa <- str_trim(temptaxa)
temptree <- unlist(strsplit(temptree," ",fixed=TRUE))
temptree <- temptree[length(temptree)]
temptree <- gsub("\\[.*?\\]","",temptree,fixed=FALSE)

temptree <- gsub(";","",temptree,fixed=TRUE)
temptree <- unlist(strsplit(temptree,","))

# putting the taxa names back into the tree
for (i in 1:(length(temptaxa))) {
for (j in 1:(length(temptree))) {
temppattern <- gsub("\\(","",temptree[j])
temppattern <- gsub("\\).*","",temppattern)
temppattern <- gsub(":.*","",temppattern)
if(temppattern==i) {
temptree[j] <- gsub(paste(i,":",sep=""),paste(temptaxa[i],":",sep=""),temptree[j],fixed=TRUE)
break
}
}
}

#extracting the internal node names
nodenames <- c(temptaxa,charsetable[1,-1])
nodenames <- nodenames[duplicated(nodenames)==FALSE]
nodenames <- nodenames[(length(temptaxa)+1):(length(nodenames))]
nodenames <- sort(as.numeric(nodenames))
recordtaxa <- matrix(NA, ncol=3,nrow=(length(nodenames)))
recordtaxa[,1] <- nodenames

#getting a list of taxa that are descendents of each internal node
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


columns <- nchar(recordtaxa[,2])-nchar(gsub(",","",recordtaxa[,2]))
recordtaxa <- recordtaxa[order(columns),]

#calculatingbranchlengths
newrecordtaxa <- NULL

for (i in 1:(dim(recordtaxa)[1]-1)) {
taxa <- unlist(strsplit(recordtaxa[i,2],","))
taxa <- gsub("(","",taxa,fixed=TRUE)
taxa <- gsub("\\).*","",taxa,fixed=FALSE)
replacetaxa <- paste(taxa,collapse=",")
replacetaxa <- paste("(",replacetaxa,")",sep="")
for (j in (i+1):dim(recordtaxa)[1]) {
recordtaxa[j,2] <- gsub(replacetaxa,recordtaxa[i,1],recordtaxa[j,2],fixed=TRUE)
}
for (j in 1:length(taxa)) {
taxadecom <- unlist(strsplit(taxa[j],":"))
temp <- cbind(recordtaxa[i,1],taxadecom[1],taxadecom[2])
newrecordtaxa <- rbind(newrecordtaxa,temp)
}
}

i <- dim(recordtaxa)[1]
taxa <- unlist(strsplit(recordtaxa[i,2],","))
taxa <- gsub("(","",taxa,fixed=TRUE)
taxa <- gsub("\\).*","",taxa,fixed=FALSE)
for (j in 1:length(taxa)) {
taxadecom <- unlist(strsplit(taxa[j],":"))
temp <- cbind(recordtaxa[i,1],taxadecom[1],taxadecom[2])
newrecordtaxa <- rbind(newrecordtaxa,temp)
}

recordtaxa <- newrecordtaxa

recordtaxa <- recordtaxa[order(as.numeric(recordtaxa[,2])),]
recordtaxa <- recordtaxa[order(as.numeric(recordtaxa[,1])),]
recordtaxa <- t(recordtaxa)
namerow <- rbind("ancest_node","daughter_node","branchlength")
recordtaxa <- cbind(namerow,recordtaxa)

taxa_plus_nodenames <- c(temptaxa,nodenames)
charsbybranchlength <- NULL

for (i in 1:(length(taxa_plus_nodenames)-1)) {
if (!(taxa_plus_nodenames[i]=="2")) {
tempbranches1 <- as.matrix(recordtaxa[1:3,(which(recordtaxa[2,]==taxa_plus_nodenames[i]))])
x <- as.numeric(tempbranches1[1,1])
while (x > 2) {
temptemp <- as.matrix(recordtaxa[1:3,suppressWarnings(which(as.numeric(recordtaxa[2,])==x))])
x <- as.numeric(temptemp[1,1])
tempbranches1 <- cbind(tempbranches1,temptemp)
}
} else {
tempbranches1 <- as.matrix(c(2, 0, 0),ncol=1)
}
for (j in (i+1):(length(taxa_plus_nodenames))) {
identicalchars <- sum(charsetable[2:(dim(charsetable)[1]),(which(charsetable[1,]==taxa_plus_nodenames[i]))]==charsetable[2:(dim(charsetable)[1]),(which(charsetable[1,]==taxa_plus_nodenames[j]))])
if (taxa_plus_nodenames[j]=="2") {
tempbranches2 <- as.matrix(c(2, 0, 0),ncol=1)
} else {
tempbranches2 <- as.matrix(recordtaxa[1:3,(which(recordtaxa[2,]==taxa_plus_nodenames[j]))])
while (!(any(tempbranches2[1,] %in% tempbranches1[1,]))) {
x <- as.numeric(tempbranches1[1,1])
temptemp <- as.matrix(recordtaxa[1:3,suppressWarnings(which(as.numeric(recordtaxa[2,])==x))])
tempbranches2 <- cbind(tempbranches1,temptemp)
}
}
sumbranchlength <- sum(as.numeric(tempbranches1[3,1:(which(tempbranches1[1,] %in% tempbranches2[1,(dim(tempbranches2)[2])]))])) + sum(as.numeric(tempbranches2[3,1:(dim(tempbranches2)[2])]))
temp <- c(taxa_plus_nodenames[i],taxa_plus_nodenames[j],identicalchars,sumbranchlength)
charsbybranchlength <- rbind(charsbybranchlength,temp)
}
}







