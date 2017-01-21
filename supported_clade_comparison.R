supported_clade_comparison <- function(path_file_to_tree1,path_file_to_tree2,cutoff) {

#path_file_to_tree1 <- "C:\\Users\\a499a400\\Dropbox\\Kaloula frogs\\Nov2015_submission\\beast\\Run 1\\16Jun2015\\annotated"
#path_file_to_tree2 <- "C:\\Users\\a499a400\\Dropbox\\Kaloula frogs\\Nov2015_submission\\beast\\Server_run_extensions\\annotated_14July2015_concatenated"
#cutoff <- 0.95

tree1 <- readLines(path_file_to_tree1)
tree2 <- readLines(path_file_to_tree2)

tree1_translate <- which(grepl("ranslate",tree1,fixed=TRUE))
tree2_translate <- which(grepl("ranslate",tree2,fixed=TRUE))

tree1_tree <- which(grepl("tree",tree1,fixed=TRUE))[2]
tree1_taxa <- unlist(strsplit(tree1[(tree1_translate+1):(tree1_tree-2)]," "))
tree1_taxa <- tree1_taxa[which(tree1_taxa!="")]
tree1_taxa <- tree1_taxa[which(tree1_taxa!="\t\t")]
tree1_taxa <- matrix(tree1_taxa,ncol=2,byrow=TRUE)
tree1_taxa[,2] <- gsub(",","",tree1_taxa[,2])

tree2_tree <- which(grepl("tree",tree2,fixed=TRUE))[2]
tree2_taxa <- unlist(strsplit(tree2[(tree2_translate+1):(tree2_tree-2)]," "))
tree2_taxa <- tree2_taxa[which(tree2_taxa!="")]
tree2_taxa <- tree2_taxa[which(tree2_taxa!="\t\t")]
tree2_taxa <- matrix(tree2_taxa,ncol=2,byrow=TRUE)
tree2_taxa[,2] <- gsub(",","",tree2_taxa[,2])

tree1 <- tree1[tree1_tree]
tree1 <- unlist(strsplit(tree1," "))[4]
tree2 <- tree2[tree2_tree]
tree2 <- unlist(strsplit(tree2," "))[4]

tree1 <- gsub("&.*?length_range=\\{.*?\\},","",tree1)
tree1 <- gsub("rate=.*?rate_range=\\{.*?\\}","",tree1)
tree1 <- gsub("rate=.*?\\]","\\]",tree1)
tree1 <- gsub("&.*?posterior","posterior",tree1)
tree1 <- unlist(strsplit(tree1,"posterior"))

tree1_output <- c("posterior","species")
  
for (i in 2:(length(tree1))) {
  posterior <- as.numeric(gsub("=","",(unlist(strsplit(tree1[i],","))[1]),fixed=TRUE))
  if(posterior>=cutoff) {
    
    # This bit needs to be tweaked - ideally you would go backwards through the brackets array until the count of "(" the count of ")" and then figure out how far this in from the beginning
    vect <- paste(tree1[1:(i-1)],collapse="")
    brackets <- unlist(strsplit(vect,"[0-9]+"))
    opens <- 0
    closes <- 0
    for (j in (length(brackets)):1) {
        opens <- opens + nchar(brackets[j]) - nchar(gsub("\\(","",brackets[j]))
        closes  <- closes + nchar(brackets[j]) - nchar(gsub("\\)","",brackets[j]))
        if (closes <= opens) {
          log_closes <- closes
        }
    }
    
    brackets <- unlist(strsplit(vect,"\\("))
    no_vect <- paste(brackets[(log_closes+1):(length(brackets))],collapse="")
    no_vect <- unlist(strsplit(no_vect,","))
    taxa <- NULL
    for (j in 1:(length(no_vect))) {
      if (length(unlist(strsplit(no_vect[j],"\\[\\]")))>1) {
        taxa <- c(taxa, unlist(strsplit(no_vect[j],"\\[\\]"))[1])
      }
    }
    temp <- c(posterior, paste(sort(tree1_taxa[(which(tree1_taxa[,1] %in% taxa)),2]),collapse=","))
    tree1_output <- rbind(tree1_output,temp)
  }
}
    
tree1_output_title <- tree1_output[1,]
tree1_output <- tree1_output[-1,]
tree1_output <- tree1_output[order(tree1_output[,2]),]
tree1_output <- rbind(tree1_output_title,tree1_output)


tree2 <- gsub("&.*?length_range=\\{.*?\\},","",tree2)
tree2 <- gsub("rate=.*?rate_range=\\{.*?\\}","",tree2)
tree2 <- gsub("rate=.*?\\]","\\]",tree2)
tree2 <- gsub("&.*?posterior","posterior",tree2)
tree2 <- unlist(strsplit(tree2,"posterior"))  
  
  
tree2_output <- c("posterior","species")
  
for (i in 2:(length(tree2))) {
  posterior <- as.numeric(gsub("=","",(unlist(strsplit(tree2[i],","))[1]),fixed=TRUE))
  if(posterior>=cutoff) {
    
    # This bit needs to be tweaked - ideally you would go backwards through the brackets array until the count of "(" the count of ")" and then figure out how far this in from the beginning
    vect <- paste(tree2[1:(i-1)],collapse="")
    brackets <- unlist(strsplit(vect,"[0-9]+"))
    opens <- 0
    closes <- 0
    for (j in (length(brackets)):1) {
        opens <- opens + nchar(brackets[j]) - nchar(gsub("\\(","",brackets[j]))
        closes  <- closes + nchar(brackets[j]) - nchar(gsub("\\)","",brackets[j]))
        if (closes <= opens) {
          log_closes <- closes
        }
    }
    
    brackets <- unlist(strsplit(vect,"\\("))
    no_vect <- paste(brackets[(log_closes+1):(length(brackets))],collapse="")
    no_vect <- unlist(strsplit(no_vect,","))
    taxa <- NULL
    for (j in 1:(length(no_vect))) {
      if (length(unlist(strsplit(no_vect[j],"\\[\\]")))>1) {
        taxa <- c(taxa, unlist(strsplit(no_vect[j],"\\[\\]"))[1])
      }
    }
    temp <- c(posterior, paste(sort(tree2_taxa[(which(tree2_taxa[,1] %in% taxa)),2]),collapse=","))
    tree2_output <- rbind(tree2_output,temp)
  }
}
    
tree2_output_title <- tree2_output[1,]
tree2_output <- tree2_output[-1,]
tree2_output <- tree2_output[order(tree2_output[,2]),]
tree2_output <- rbind(tree2_output_title,tree2_output)  
  
  
print(paste(sum(tree1_output[,2] %in% tree2_output[,2])," clades out of ",dim(tree1_output)[1]," highly supported clades in tree 1 are also found in tree 2.",sep=""))
print(paste(sum(tree2_output[,2] %in% tree1_output[,2])," clades out of ",dim(tree2_output)[1]," highly supported clades in tree 2 are also found in tree 1.",sep="")) 

print("The following clades were highly supported in tree 1 but not highly supported or not found in tree 2")
print(rbind(tree1_output_title,tree1_output[(which(!(tree1_output[,2] %in% tree2_output[,2]))),]))
  
print("The following clades were highly supported in tree 2 but not highly supported or not found in tree 1")
print(rbind(tree2_output_title,tree2_output[(which(!(tree2_output[,2] %in% tree1_output[,2]))),]))

}
