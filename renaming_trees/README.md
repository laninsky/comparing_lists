Change working directory to whatever is appropriate:
```
setwd("C:\\Users\\a499a400\\Downloads")
namelist <- as.matrix(read.table("namelist.txt"))
treefilename <- "cornufer.nex"

temp <- readLines(treefilename)
lentemp <- length(temp)
lennamelist <- dim(namelist)[1]
outtemp <- matrix(NA,nrow=lentemp,ncol=1)
