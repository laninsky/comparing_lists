# renaming_things_inside_nexus_file

Things you need:

-- A tab-delimited txt document with the names you want in the left-hand column (namelist.txt), and the names you want replaced in the second (an example given with the scripts above). These names can't have any white space in them (either the names you have, or the names you want). If the names you have already have spaces in them, replace them in the .nex file/namelist.txt with sed -i.

-- Your nexus files you want to have the names replaced within (note: these should be sequence files only! For tree files represented as a single nexus line, see the folder "renaming_trees").

-- For the nexus files, namelist.txt, and the script rename.R to be in the same directory. Once they are, you can cycle through multiple nexus files by copying and pasting the script in rename.sh. The R script assumes you have previously installed the R package 'stringr'.

This pipeline wouldn't be possible without:

R: R Core Team. 2015. R: A language and environment for statistical computing. URL http://www.R-project.org/. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/

Stringr:  Hadley Wickham (2012). stringr: Make it easier to work with strings..
  R package version 0.6.2. http://CRAN.R-project.org/package=stringr (for up-to-date citation information run citation("stringr" in R)

