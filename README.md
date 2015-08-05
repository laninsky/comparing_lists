# running_SVDquartets
Running SVDquartets on your *.nex file (similar to what SNAPP eats)

#Step 1
Install PAUP (http://people.sc.fsu.edu/~dswofford/paup_test/) and then run it from the directory with your *.nex files:
```
/scratch/a499a400/bin/paup4a146_centos64
```

#Step 2
Execute your data e.g.
```
execute SNAPP_60.nex;
```
If there are any issues make sure none of your sample names have hyphens in them.

#Step 3
```
SVDQuartets evalQuartets=all showScores=yes qfile=SNAPP_60.qfile;                 
end;
quit;
```
Make sure to change the tree file to whatever you would rather!

#Step 4 -- THIS IS AS FAR AS I HAVE GOT
Get a copy of Quartet MaxCut (http://research.haifa.ac.il/~ssagi/software/QMCN.tar.gz)
````
/scratch/a499a400/bin/QMCN/find-cut-Linux-64 qrtt=SNAPP_60.qfile
/scratch/a499400/bin/QMCN/quartet-agreement-Linux-64 tree=MXCNtree.dat qrtt=SNAPP_60.qfile 
/scratch/a499a400/bin/QMCN/genTreeAndQuartets-Linux-64 42 111930 55.9 0


