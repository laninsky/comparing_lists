# converting_concat_phylip_to_fasta_SNP
Originally at: https://github.com/laninsky/converting_concat_phylip_to_fasta_SNP

```
setwd("C:\\Users\\a499a400\\Dropbox\\Kaloula frogs")
intable <- readLines("froggies-25pctincomplete-filled-MAFFT-gblocks.nexus.phylip")
notaxa <- as.numeric(unlist(strsplit(intable[1],"\\s+"))[1])
noseq <- as.numeric(unlist(strsplit(intable[1],"\\s+"))[2])

seqmatrix <- matrix(NA,ncol=noseq,nrow=notaxa)

for (i in 1:notaxa) {
concatseq <- unlist(strsplit(intable[i+1],"\\s+"))[2]
seqmatrix[i,] <- unlist(strsplit(concatseq,""))
}

colstokeep <- matrix(NA,ncol=noseq,nrow=1)
for (i in 1:noseq) {
colstokeep[1,i] <- sum(unique(seqmatrix[,i])!="?")
}

colstokeepindex <- which(colstokeep[1,]>=2)

seqmatrix <- seqmatrix[,colstokeepindex]

finaloutput <-  matrix(NA,ncol=1,nrow=(notaxa*2))
nums <- (seq(1,(notaxa*2),2))

for (i in 1:length(nums)) {
finaloutput[nums[i],1] <- paste(">",(unlist(strsplit(intable[(i+1)],"\\s+"))[1]),sep="")
finaloutput[(nums[i]+1),1] <- paste(seqmatrix[i,],collapse="")
}

write.table(finaloutput, "varsitesonly.fas",quote=FALSE, col.names=FALSE,row.names=FALSE)
```

### Version history
This script was written for the following manuscript:
Alexander, A.M., Su, Y.C., Oliveros, C.H., Olson, K.V., Travers, S.L. and Brown, R.M., 2017. Genomic data reveals potential for hybridization, introgression, and incomplete lineage sorting to confound phylogenetic relationships in an adaptive radiation of narrow‐mouth frogs. Evolution, 71(2), pp.475-488. 

I am no longer actively maintaining this repository, but will respond to issues.

This script wouldn't be possible without:  
R: R Core Team. 2015. R: A language and environment for statistical computing. URL http://www.R-project.org/. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/
