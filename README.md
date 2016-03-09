#Mapping_SNPs_to_trees
You want to assess how well your SNPs "fit" a given tree

So the way I have gone about this is to use MESQUITE. I converted a sequential Phylip file to a nexus file (with invariant sites excluded) with the code below:
```
library(stringr)
setwd("working_dir") # e.g. setwd("C:/Users/a499a400/Dropbox/Kaloula frogs")
input <- readLines("filename") # e.g. input <- readLines("froggies-25pctincomplete-filled-MAFFT-gblocks.nexus.phylip")
notaxa <- as.numeric(unlist(strsplit(input[1], " "))[1])
nosites <-  as.numeric(unlist(strsplit(input[1], " "))[2])
inputmatrix <- matrix(ncol=(nosites+1),nrow=notaxa)

for (i in 1:notaxa) {
temp <- unlist(strsplit(input[i+1]," "))
inputmatrix[i,1] <- temp[1]
inputmatrix[i,2:(dim(inputmatrix)[2])] <- unlist(strsplit(temp[2],""))
}

states <- NULL
for (i in 1:nosites) {
states[i] <- sum(unique(inputmatrix[,(i+1)])!="?")
}

outputmatrix <- inputmatrix[,(which(states[1:nosites]>1)+1)]
nosites <- dim(outputmatrix)[2]

tempmatrix <- matrix(ncol=2,nrow=notaxa)
tempmatrix[,1] <- inputmatrix[,1]
for (i in 1:notaxa) {
tempmatrix[i,2] <- paste(outputmatrix[i,1:(dim(outputmatrix)[2])],collapse="")
}
outputmatrix <- matrix(ncol=1,nrow=notaxa)
for (i in 1:notaxa) {
outputmatrix[i,1] <- paste(tempmatrix[i,],collapse=" ")
}
outputmatrix <- rbind("matrix",outputmatrix)
outputmatrix <- rbind("  format datatype=dna missing=? gap=-;",outputmatrix)
temp <- paste("  dimensions ntax=",notaxa," nchar=",nosites,";",sep="")
outputmatrix <- rbind(temp,outputmatrix)
outputmatrix <- rbind("begin data;",outputmatrix)
outputmatrix <- rbind("#NEXUS",outputmatrix)
outputmatrix <- rbind(outputmatrix,";")
outputmatrix <- rbind(outputmatrix,"end;")

write.table(outputmatrix, "output.nexus",quote=FALSE, col.names=FALSE,row.names=FALSE)
```

I imported this as the matrix (MESQUITE is funny about taxa names, so make sure you have alphabet characters as well as sample numbers), and imported the associated tree (with identical taxa names). In the tree window, I then went to Analysis:Tree > Trace all characters. This spits out a text window with the inferred character at each internal node. In this text window, I then selected "show terminal taxa" to make sure we got everything in there, copied it out, and pasted it into word.

In word, I deleted the top rows describing the dataset (keeping everything from 'Char.\Node" onwards). I went through word, because I wanted to keep the tab delimited structure, because some of the characters are linked e.g. A G vs A, and therefore have a whitespace between them. After deleting these lines, I saved this as a txt file. Mesquite numbers the internal nodes in order that they appear in the tree file.

The next step is using R to work out state changes in the SNPs between nodes. To do this, copy/source the R script (SNPs_to_clades.R) in this repository. Invoke it by: 
```
SNPs_to_clades(working_dir,charset,tree) 
```
Where your charset is the output from MESQUITE you edited in word to remove the top lines, and your tree is the same tree as you imported into MESQUITE e.g. 
```
SNPs_to_clades("C:/Users/a499a400/Dropbox/chan","Char.txt","ExaBayes_ConsensusExtendedMajorityRuleNexus.contree.tre")
```
#Output
You should a two tab-delimited file as output (state_changes_along_branches.txt). This file has three rows at the top - the first row has the ancestral nodes, and the second row has the daughter nodes (in this file, each of the daughters has its own column, rather than being comma separated), and the third row has the branchlengths between these. In the left hand column of the subsequent rows are the character names (e.g. SNPs). For each column of ancestor > daughter, the number of inferred state changes along the branch, divided by branchlength, are given for each SNP.

#MSDS
You are then probably going to want to visualize branches in the tree which seem to have more or less state changes given their length. These might be longer branch lengths (where potentially multiple substitutions have masked each other), or potentially areas of the genome affected by introgression etc.
```
setwd("C:/Users/a499a400/Dropbox/chan")
recordtaxa <- as.matrix(read.table("state_changes_along_branches.txt",sep="\t"))

d <- dist(t(recordtaxa[3:(dim(recordtaxa)[1]),2:(dim(recordtaxa)[2])]))
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric	MDS",	type="n")
text(x, y, cex=.7)
```
In the resulting plot, you might see some numbers as outliers. These are your problem branches if so. Create an array with them:
```
problem <- c(119,59,80)
problem <- problem+1
recordtaxa[1:2,problem]
```
You probably want to then exclude these outliers and run the plot again, just in case they were masking further outliers.
```
recordtaxa <- recordtaxa[,-problem]
```
Head back up and run the d <- dist function again etc.

#Correlations in changes between branches

```
setwd("C:/Users/a499a400/Dropbox/chan")
recordtaxa <- as.matrix(read.table("state_changes_along_branches.txt",sep="\t"))

rsquares <- NULL

for (i in 2:(dim(recordtaxa)[2]-1)) {
for (j in (i+1):(dim(recordtaxa)[2])) {
lmsum <- lm(recordtaxa[3:(dim(recordtaxa)[1]),i] ~ recordtaxa[3:(dim(recordtaxa)[1]),j])
temp <- c(i, j, summary(lmsum)$r.squared)
rsquares <- rbind(rsquares, temp)
}
}

rsquares <- rsquares[order(rsquares[,3]),]
```
For the pairs of branches that show a strong enough correlation in which SNPs have changes along the, for you to be interested in them, you can query them by:
```
recordtaxa[1:2,(rsquares[rowofinterest,1])]
recordtaxa[1:2,(rsquares[rowofinterest,2])]
```
e.g. I was interested in the last row of rsquares:
```
108  119 6.582298e-01
recordtaxa[1:2,(rsquares[(dim(rsquares)[1]),1])]
[1] "1.070000e+02" "1.190000e+02"
recordtaxa[1:2,(rsquares[(dim(rsquares)[1]),2])]
[1] "119.00000000" "121.00000000"
```
i.e. in this dataset, the branch from internal nodes 107 > 119 appears to show a lot of correlated changes in specific SNPs with the branch from internal nodes 119 > 121. Interestingly, in this dataset, there is a known hybrid that shares ancestral node 119 with the branch that leads to 121.


