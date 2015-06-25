# renaming_things_inside_nexus_file

Things you need:

-- A tab-delimited txt document with the names you want in the left-hand column (namelist.txt), and the names you want replaced in the second (an example given with the scripts above). These names can't have any white space in them (either the names you have, or the names you want). If the names you have already have spaces in them, replace them in the .nex file/namelist.txt with sed -i.

-- Your nexus files you want to have the names replaced within.

-- For the nexus files, namelist.txt, and the script rename.R to be in the same directory. Once they are, you can cycle through multiple nexus files by copying and pasting the script in rename.sh. The R script assumes you have previously installed the R package 'stringr'.
