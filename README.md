# running_SVDquartets
Running SVDquartets on your *.nex file (similar to what SNAPP eats e.g. sites should be unlinked)

#Step 1
Install PAUP (http://people.sc.fsu.edu/~dswofford/paup_test/) and then run it from the directory with your *.nex files:

#Step 2
Modify the PBS file below to execute it for your data e.g.

```
/scratch/a499a400/bin/paup4a150_centos64 
execute SNAPP_phased.nex; SVDQuartets evalQuartets=all partition=frogspecies speciesTree=no showScores=yes qfile=SNAPP_60.qfile savetrees; SVDQuartets evalQuartets=all speciesTree=no showScores=yes bootstrap nreps=500 nthreads=6; savetrees file=SNAPP_boot.qfm;
```
The 50% majority consensus bootstrap tree will be printed to SNAPP_boot.qfm, actual tree will be SNAPP.qfm. The bootstraps are rendered as 'Branch Times' in Fig tree.

The svdq ? command gives you a list of the available options for use with SVDQuartets if you run it after launching paup:
```
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
```

For more information on this, check out:
http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2015/01/SVDquartets_tutorial2015.pdf

http://www.stat.osu.edu/~lkubatko/SVDquartets_tutorial2015.html

If you don't evaluate all of the quartets, bump up the number of quartets sampled

Make sure to change the qfile, bootstrap treeFile and tree file to whatever you would rather call them. This method is going to output the QFM tree inference method.

