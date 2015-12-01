# To get a list of overlapping UCE loci between two different datasets
Get a list of loci from the incomplete folders for each of your taxa e,g,:

```
ls *.fasta > kaloula_75_uces.txt
ls *.nexus > sanguirana_75_uces.txt
```

Then, in R:

```
working_dir <- "C:/Users/a499a400/Dropbox/Kaloula frogs/uce_overlaps/"
file1 <- "kaloula_taxa_list.txt"
file2 <- "sang_taxa_list.txt"

setwd(working_dir)
tempfile1 <- read.table(file1,sep=".")
tempfile2 <- read.table(file2,sep=".")

loci_found_in_both <- tempfile1[(tempfile1[,1] %in% tempfile2[,1]),1]

write.table(loci_found_in_both,"loci_found_in_both.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
```

#To get a list of overlappng UCE loci between two different datasets, and the number of pis in each species (following https://github.com/laninsky/No_of_pis_per_locus)

In R:
```
working_dir <- "C:/Users/a499a400/Dropbox/Kaloula frogs/uce_overlaps/"
file1 <- "list_of_pis_by_locus_kaloula.txt"
file2 <- "list_of_pis_by_locus_sanguirana.txt"

setwd(working_dir)
tempfile1 <- as.matrix(read.table(file1,sep=" "))
tempfile2 <- as.matrix(read.table(file2,sep=" "))

tempfile1big <- matrix("", ncol=4,nrow=(dim(tempfile1)[1]))
tempfile2big <- matrix("", ncol=4,nrow=(dim(tempfile2)[1]))

for (j in 1:dim(tempfile1)[1]) {
tempfile1big[j,1] <- unlist(strsplit(tempfile1[j,1],".", fixed=TRUE))[1]
}

for (j in 1:dim(tempfile2)[1]) {
tempfile2big[j,1] <- unlist(strsplit(tempfile2[j,1],".", fixed=TRUE))[1]
}

tempfile1big[,2] <- tempfile1[,2]
tempfile1big[,3] <- tempfile1[,3]
suppressWarnings(tempfile1big[,4] <- as.numeric(tempfile1[,2])/as.numeric(tempfile1[,3]))
tempfile1big[1,4] <- "pis/length"

tempfile2big[,2] <- tempfile2[,2]
tempfile2big[,3] <- tempfile2[,3]
suppressWarnings(tempfile2big[,4] <- as.numeric(tempfile2[,2])/as.numeric(tempfile2[,3]))
tempfile2big[1,4] <- "pis/length"

loci_found_in_both <- tempfile1[(tempfile1[,1] %in% tempfile2[,1]),1]

header <- c("locus","sp1_pis","sp1_length","sp1_pis/length","sp2_pis","sp2_length","sp2_pis/length")

output <- matrix("", ncol=7,nrow=(length(loci_found_in_both)-1))

output[,1] <- loci_found_in_both[2:(length(loci_found_in_both))]

x <- 1
for (i in 2:(dim(tempfile1big)[1])) {
tempname <- paste(tempfile1big[i,1],".nexus",sep="")
if (tempname==output[x,1]) {
output[x,2:4] <- tempfile1big[i,2:4]
x <- x+1
}
}

x <- 1
for (i in 2:(dim(tempfile2big)[1])) {
tempname <- paste(tempfile2big[i,1],".nexus",sep="")
if (tempname==output[x,1]) {
output[x,5:7] <- tempfile2big[i,2:4]
x <- x+1
}
}

output <- rbind(header,output)

write.table(output,"loci_found_in_both_with_pis.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
```



