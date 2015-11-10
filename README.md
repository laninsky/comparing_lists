# running_SVDquartets
Running SVDquartets on your *.nex file (similar to what SNAPP eats)

#Step 1
Install PAUP (http://people.sc.fsu.edu/~dswofford/paup_test/) and then run it from the directory with your *.nex files:

#Step 2
Modify the PBS file below to execute it for your data e.g.

```
/scratch/a499a400/bin/paup4a146_centos64 
execute SNAPP_60.nex; SVDQuartets evalQuartets=all showScores=yes qfile=SNAPP_60.qfile bootstrap treeFile=SNAPP_60_boot.qfm; savetrees file=SNAPP_60.qfm;

gettrees file=SNAPP_60_boot.qfm; contree all/strict=no majrule=yes usetreewts=yes treefile=SNAPP_60_con.tre;

execute SNAPP_80.nex; SVDQuartets evalQuartets=all showScores=yes qfile=SNAPP_80.qfile bootstrap treeFile=SNAPP_80_boot.qfm; savetrees file=SNAPP_80.qfm; gettrees file=SNAPP_80_boot.qfm; 

gettrees file=SNAPP_80_boot.qfm; contree all/strict=no majrule=yes usetreewts=yes treefile=SNAPP_80_con.tre; quit;

```

The svdq ? command gives you a list of the available options for use with SVDQuartets if you run it after launching paup:
Keyword ------- Option type --------------------- Current setting ----------

evalQuartets    all|random                        random

nquartets       <real-value>                      100000

preferAllQ      no|yes                            yes

speciesTree     no|yes                            no

partition       <taxpartition-name>               (none)

treeInf         QFM|curTrees|none                 QFM

seed            <integer-value>                   0

bootstrap       no|yes                            no

nreps           <integer-value>                   100

nthreads        ncpus|<number-of-threads>         2

mrpFile         <species-outfile-name>            (none)

qfile           <quartets-outfile-name>           (none)

qformat         qmc|qfm                           qmc

replace         no|yes                           *no

showScores      no|yes                            no

showSV          no|yes                            no

treeFile        <filename-for-bootstrap-treefile> (none)

treemodel       mscoalescent|shared               mscoalescent
                                                 *Option is nonpersistent


For more information on this, check out:
http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2015/01/SVDquartets_tutorial2015.pdf

http://www.stat.osu.edu/~lkubatko/SVDquartets_tutorial2015.html

If you don't evaluate all of the quartets, bump up the number of quartets sampled

Make sure to change the qfile, bootstrap treeFile and tree file to whatever you would rather call them



#Step 4 -- THIS IS AS FAR AS I HAVE GOT
Get a copy of Quartet MaxCut (http://research.haifa.ac.il/~ssagi/software/QMCN.tar.gz)
````
/scratch/a499a400/bin/QMCN/find-cut-Linux-64 qrtt=SNAPP_60.qfile
/scratch/a499400/bin/QMCN/quartet-agreement-Linux-64 tree=MXCNtree.dat qrtt=SNAPP_60.qfile 
/scratch/a499a400/bin/QMCN/genTreeAndQuartets-Linux-64 42 111930 55.9 0



