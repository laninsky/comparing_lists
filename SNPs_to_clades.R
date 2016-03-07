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

temptree

