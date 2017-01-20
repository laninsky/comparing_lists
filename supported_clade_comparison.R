supported_clade_comparison(path_file_to_tree1,path_file_to_tree2,cutoff) {

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
tree2 <- gsub("&.*?length_range=\\{.*?\\},","",tree2)
tree2 <- gsub("rate=.*?rate_range=\\{.*?\\}","",tree2)
