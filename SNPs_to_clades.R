SNPs_to_clades <- function(working_dir,charset,tree) {

#loading libraries and loading in files
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

# making the tree and taxa variables prettier
temptaxa <- temptaxa[-1,1]
subsample <- seq(3,(length(temptaxa)*3),3)
temptaxa <- unlist(strsplit(temptaxa,"\t",fixed=TRUE))[subsample]
temptree <- unlist(strsplit(temptree," ",fixed=TRUE))
temptree <- temptree[length(temptree)]
temptree <- gsub("\\[.*?\\]","",temptree,fixed=FALSE)
branchlengths <- unlist(strsplit(temptree,":",fixed=TRUE))
branchlengths <- branchlengths[-1]
branchlengths <- gsub(").*","",branchlengths,fixed=FALSE)
branchlengths <- gsub("\\(.*","",branchlengths,fixed=FALSE)
branchlengths <- gsub(",.*","",branchlengths,fixed=FALSE)

temptree <- gsub("\\:[0-9]+.[0-9]+[e-]*[0-9]*","",temptree,fixed=FALSE)
temptree <- gsub(";","",temptree,fixed=TRUE)
temptree <- unlist(strsplit(temptree,","))

# putting the taxa names back into the tree
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

#extracting the internal node names
nodenames <- c(temptaxa,charsetable[1,-1])
nodenames <- nodenames[duplicated(nodenames)==FALSE]
nodenames <- nodenames[(length(temptaxa)+1):(length(nodenames))]
nodenames <- sort(as.numeric(nodenames))
recordtaxa <- matrix(NA, ncol=2,nrow=(length(nodenames)))
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

#substituting internal node names for groups of taxa (so  
recordtaxa[1:dim(recordtaxa)[1],2] <- gsub("(","",recordtaxa[1:dim(recordtaxa)[1],2],fixed=TRUE)
recordtaxa[1:dim(recordtaxa)[1],2] <- gsub(")","",recordtaxa[1:dim(recordtaxa)[1],2],fixed=TRUE)
recordtaxa <- recordtaxa[order(nchar(recordtaxa[,2])),]
for (i in 1:((dim(recordtaxa)[1])-1)) {
for (j in (i+1):(dim(recordtaxa)[1])) {
recordtaxa[j,2] <- gsub(recordtaxa[i,2],recordtaxa[i,1],recordtaxa[j,2],fixed=TRUE)
}
}

recordtaxa <- recordtaxa[order(as.numeric(recordtaxa[,1])),]

#adding a column for branchlengths
rowadd <- matrix(NA,nrow=dim(recordtaxa)[1],ncol=1)
recordtaxa <- cbind(recordtaxa,rowadd)

x <- 1
for (i in 1:dim(recordtaxa)[1]) {
numbranches <- (nchar(recordtaxa[i,2]) - nchar(gsub(",","",recordtaxa[i,2],fixed=TRUE)))
recordtaxa[i,3] <- paste(branchlengths[x:(x+numbranches)],collapse=",")
x <- numbranches + x + 1
}

write.table(recordtaxa,"list_of_clades_in_tree.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

# giving each daughter node its own row
for (i in 1:(dim(recordtaxa)[1])) {
newsecond <- unlist(strsplit(recordtaxa[i,2],","))
newthird <- unlist(strsplit(recordtaxa[i,3],","))
recordtaxa[i,2] <- newsecond[1]
recordtaxa[i,3] <- newthird[1]
for (j in 2:(length(newsecond))) {
temprow <- cbind(recordtaxa[i,1],newsecond[j],newthird[j])
recordtaxa <- rbind(recordtaxa,temprow)
}
}

recordtaxa <- recordtaxa[order(as.numeric(recordtaxa[,2])),]
recordtaxa <- recordtaxa[order(as.numeric(recordtaxa[,1])),]
recordtaxa <- t(recordtaxa)
namerow <- rbind("ancest_node","daughter_node","branchlength")
recordtaxa <- cbind(namerow,recordtaxa)

newmatrix <- matrix(NA,nrow=(dim(charsetable)[1]-1),ncol=(dim(recordtaxa)[2]))
recordtaxa <- rbind(recordtaxa,newmatrix)
recordtaxa[4:(dim(recordtaxa)[1]),1] <- charsetable[2:(dim(charsetable)[1]),1]

for (i in 4:(dim(recordtaxa)[1])) {
for (j in 2:(dim(recordtaxa)[2])) {
first <- which((charsetable[1,]==recordtaxa[1,j])==TRUE)
second <- which((charsetable[1,]==recordtaxa[2,j])==TRUE)

if ((charsetable[(i-2),first]=="N" | charsetable[(i-2),second]=="N")==TRUE) {
break
}
if(charsetable[(i-2),first]==charsetable[(i-2),second]) {
recordtaxa[i,j] <- 0
} else {
recordtaxa[i,j] <- 1/(as.numeric(recordtaxa[3,j]))
}
}
}

write.table(recordtaxa,"state_changes_along_branches.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

}
