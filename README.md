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
