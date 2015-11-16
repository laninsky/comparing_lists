setwd("C:/Users/a499a400/Dropbox/Kaloula frogs/Genbank_submission")

temp <- as.matrix(read.table("Raw_data.txt",sep="\t"))

names <- temp[1,(seq(8,77,3))]

featuretable <- NULL
sequence <- NULL
tempfeature <- NULL

seqs <- (seq(8,77,3))

for (t in 1:24) {
i <- seqs[t]
x <- 1
seqidentifier <- paste(">",names[t],"_",x,"[organism=Kaloula baleata] [molecule=DNA] [location=mitochondrion]",sep="")
sequence <- rbind(sequence,seqidentifier)
featuretable <- rbind(featuretable, (paste(">Features ",names[t],"_",x,sep="")))
record <- 1
tempsequence <- NULL

for (j in 2:(dim(temp)[1])) {

if(temp[j,(i-2)]=="N") {
if(!(temp[(j-1),(i-2)]=="N")) {
x <- x + 1
sequence <- rbind(sequence,tempsequence)
featuretable <- rbind(featuretable, (paste(">Features ",names[t],"_",x,sep="")))
seqidentifier <- paste(">",names[t],"_",x,"[organism=Kaloula baleata] [molecule=DNA] [location=mitochondrion]",sep="")
sequence <- rbind(sequence,seqidentifier)
record <- 2
tempsequence <- NULL
}
} else {
if(!(temp[j,(i-2)]=="-")) {
if(record==2) {
tempfeature <- paste(j,temp[j,4])
featuretable <- rbind(featuretable,tempfeature)
record <- 1
}
if(!(temp[(j-1),4]==temp[j,4])) {
tempfeature <- paste(j,temp[j,4])
featuretable <- rbind(featuretable,tempfeature)
}
tempsequence <- paste(tempsequence,temp[j,i],sep="")
}
}
}

sequence <- rbind(sequence,tempsequence)
}

write.table(sequence,"sequence.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

write.table(featuretable,"featuretable.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
