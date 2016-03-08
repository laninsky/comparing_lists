#Mapping_SNPs_to_trees
You don't want to make a ton of gene trees, but you do want to assess how well your SNPs "fit" a given tree

So the way I have gone about this is to use MESQUITE. I converted a Phylip SNP file to a nexus file (through MEGA) and imported this as the matrix (MESQUITE is funny about taxa names, so I inserted a t in front of all the sample numbers when it was a MEGA file, before converting to a nexus file), and imported the associated tree (also putting a t in front of each taxa name). In the tree window, I then went to Analysis:Tree > Trace all characters. This spits out a text window with the inferred character at each internal node. In this text window, I then selected "show terminal taxa" to make sure we got everything in there, copied it out, and pasted it into word.

In word, I deleted the top rows describing the dataset (keeping everything from 'Char.\Node" onwards). I went through word, because I wanted to keep the tab delimited structure, because some of the characters are linked e.g. A G vs A, and therefore have a whitespace between them. After deleting these lines, I saved this as a txt file. Mesquite numbers the internal nodes in order that they appear in the tree file.

The next step is using R to work out state changes in the SNPs between nodes. To do this, copy/source the R script (SNPs_to_clades.R) in this repository. Invoke it by: 

SNPs_to_clades(working_dir,charset,tree) 

Where your charset is the output from MESQUITE you edited in word to remove the top lines, and your tree is the same tree as you imported into MESQUITE

e.g. SNPs_to_clades("C:/Users/a499a400/Dropbox/chan","Char.txt","ExaBayes_ConsensusExtendedMajorityRuleNexus.contree.tre")

#Output
You should get two tab-delimited files as output: the first (list_of_clades_in_tree.txt), lists ancestral daughter nodes in the left column, and daughter nodes of these in the middle column (separated by commas), and the branchlengths between the ancestor and daughter nodes (also tab delimited) in the right column. This file is just produced for your reference.

The second file is probably the one you will be more interested in (state_changes_along_branches.txt). This file has three rows at the top - the first row has the ancestral node, and the second row has the daughter node (in this file, each of the daughters has its own column, rather than being comma separated), and the third row has the branchlengths between these. In the left hand column of the subsequent rows are the character names (e.g. SNPs). For each column of ancestor > daughter, the number of inferred state changes along the branch, divided by branchlength, are given for each SNP.

#MSDS
You are then probably going to want to visualize branches in the tree which seem to have more or less state changes given their length. These might be longer branch lengths (where potentially multiple substitutions have masked each other), or potentially areas of the genome affected by introgression etc.
```
setwd("C:/Users/a499a400/Dropbox/chan")
recordtaxa <- as.matrix(read.table("state_changes_along_branches.txt",sep="\t"))

d <- dist(t(recordtaxa[3:(dim(recordtaxa)[1]),2:(dim(recordtaxa)[2])]))
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric	MDS",	type="n")
text(x, y, cex=.7)
```
In the resulting plot, you might see some numbers as outliers. These are your problem branches if so. Create an array with them:
```
problem <- c(112,121)
problem <- problem+1
recordtaxa[1:2,problem]
```
You probably want to then exclude these outliers and run the plot again, just in case they were masking further outliers.
```
recordtaxa <- recordtaxa[,-problem]
```
Head back up and run the d <- dist function again etc.
