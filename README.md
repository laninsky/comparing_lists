# converting_concat_phylip_to_fasta_SNP

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
