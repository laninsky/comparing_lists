for i in `ls *.nex`;
do rm -rf temp*;
mv $i temp;
Rscript rename.R;
mv tempout $i;
done;
