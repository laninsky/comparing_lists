To see what might lead to correlations in the SNPs changing between branches, I simulated a DNA dataset using fastsimcoal. I used the same number of loci and average length as the empirical study (scenario.def and scenario.tpl in this folder), and based the times for merging off of the concatenated tree branchlengths. After executing with the command:
```
fsc25221.exe -t scenario.tpl -n 1 -f scenario.def -q > log.log
```
I then converted the arlequin file to a nexus file to be used in MESQUITE by (modify the setwd, arl and treeset locations - also this is assuming you have just one sample from each population, and that the samples were defined in the same order in the fastsimcoal file as in your treefile):
```
setwd("C:/Users/a499a400/Dropbox/Kaloula frogs/fastimcoal/scenario")
arl <- readLines("scenario_1_1.arp")
treeset <- readLines("C:/Users/a499a400/Dropbox/Kaloula frogs/annotated_14July2015_concatenated.tre")

library(stringr)

nosites <- arl[grepl("#Reporting status of a maximum of ",arl)]

if (length(nosites)==0) {
nosites <- arl[grepl("#Total number of polymorphic sites:",arl)]
nosites <- unlist(strsplit(nosites,"\\s+"))
nosites <- nosites[(length(nosites))]
} else {
nosites <- unlist(strsplit(nosites,"\\s+"))
nosites <- nosites[(length(nosites)-1)]
nosites2 <- arl[grepl("#Total number of polymorphic sites:",arl)]
nosites2 <- unlist(strsplit(nosites2,"\\s+"))
nosites2 <- nosites2[(length(nosites2))]
if (nosites2 < nosites) {
nosites <- nosites2
}
}

sequence <- arl[(which(grepl("SampleData=",arl)))+1]
for (i in 1:(length(sequence))) {
temp <- unlist(strsplit(sequence[i],"\\s+"))
sequence[i] <- temp[length(temp)]
}

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
temptaxa <- temptaxa[-1,1]
temptaxa <- gsub("\t"," ",temptaxa,fixed=TRUE)
temptaxa <- str_trim(temptaxa)

notaxa <- length(temptaxa)
tempmatrix <- matrix(ncol=2,nrow=notaxa)
tempmatrix[,1] <- temptaxa
tempmatrix[,2] <- sequence

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

write.table(outputmatrix, "simulated.nexus",quote=FALSE, col.names=FALSE,row.names=FALSE)
```
