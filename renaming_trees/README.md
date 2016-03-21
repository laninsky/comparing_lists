Change working directory to whatever is appropriate:
```
setwd("C:\\Users\\a499a400\\Downloads")
namelist <- as.matrix(read.table("namelist.txt"))
treefilename <- "cornufer.nex"

temp <- readLines(treefilename)
lentemp <- length(temp)
lennamelist <- dim(namelist)[1]
outtemp <- matrix(NA,nrow=lentemp,ncol=1)

for (i in 1:lentemp) {
if (!(grepl(namelist[1,2],temp[i],fixed=TRUE))) {
outtemp[i] <- temp[i] } else {
for (j in 1:lennamelist) {
temp[i] <- gsub(namelist[j,2],namelist[j,1],temp[i],fixed=TRUE)
}
outtemp[i] <- temp[i]
}
}


write.table(outtemp, "tempout",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
