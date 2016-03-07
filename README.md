# mapping_SNPs_to_trees
You don't want to make a ton of gene trees, but you do want to assess how well your SNPs "fit" a given tree

So the way I have gone about this is to use MESQUITE. I converted a Phylip SNP file to a nexus file (through MEGA) and imported this as the matrix (MESQUITE is funny about taxa names, so I inserted a t in front of all the sample numbers when it was a MEGA file, before converting to a nexus file), and did this with the associated tree (also putting a t in front of each taxa name). In the tree window, I then went to Analysis:Tree > Trace all characters. This spits out a text window with the inferred character at each internal node. In this text window, I then selected "show terminal taxa" to make sure we got everything in there, copied it out, and pasted it into word.

In word, I deleted the top rows describing the dataset (keeping everything from 'Char.\Node" onwards). I went through word, because I wanted to keep the tab delimited structure, because some of the characters are linked e.g. A G vs A. After deleting these lines, I saved this as a txt file. Mesquite numbers the internal nodes in order that they appear in the tree file.
