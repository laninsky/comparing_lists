structure_into_treemix <- function(working_dir,file_name) {
#example_use: structure_into_treemix("C:/Users/a499a400/Downloads","sangui_structure_popmap.txt")

setwd(working_dir)
structure <- as.matrix(read.table(file_name))
popnames <- unique(structure[,2])
noloci <- dim(structure)[2]-2
nopops <- length(popnames)
output <- matrix(NA,ncol=nopops,nrow=(noloci+1))
output[1,] <- popnames

for (i in 3:(dim(structure)[2])) {
states <- unique(structure[,i])[(unique(structure[,i])!=0)==TRUE]

if (is.null(states)) {
output[(i-1),] <- "0,0"
}

if(length(states)>2) {
print(i-2)
stop("SNP printed above has more than 2 states. Filter your file and try again")
}

if(length(states)==2) {
for (j in 1:nopops) {
state1count <- sum(structure[,2]==popnames[j] & structure[,i]==states[1])
state2count <- sum(structure[,2]==popnames[j] & structure[,i]==states[2])
tempoutput <- paste(state1count,",",state2count,sep="")
output[(i-1),j] <- tempoutput
}
}

if(length(states)==1) {
for (j in 1:nopops) {
state1count <- sum(structure[,2]==popnames[j] & structure[,i]==states[1])
tempoutput <- paste(state1count,",",0,sep="")
output[(i-1),j] <- tempoutput
}
}
}

write.table(output,"treemix_output.txt", sep=" ",quote=FALSE, row.names=FALSE,col.names=FALSE)

}
