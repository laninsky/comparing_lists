# structure_into_treemix
Have you got a structure file that you would like to convert to treemix? Step right up!

#What do I need?
You need a structure file, which has two rows of genotype SNP data per individual. This program expects missing data to be encoded as 0, and SNPs A-T as 1-4. It expects NO header file, and for the SNPs to already be filtered down to unlinked SNPs (e.g. a single SNP per locus). It also expects two initial columns, the first with your individual_IDs, and the second with your pop_IDs.

What if I have linked SNPs that I need to filter out?

You can follow the code at: https://github.com/laninsky/ambigoos_into_structure#what-if-you-want-to-tweak-the-individuals-in-the-filechange-completeness-of-datasetsnp-selection-criteria

What if I don't have pop_IDs?

You can add the pop_IDs using the code at: https://github.com/laninsky/ambigoos_into_structure#what-if-you-wanted-to-add-a-population-identifier-column-after-the-individual-column

What if the way I encode SNPs is a little different?

You can crib the code at: https://github.com/laninsky/creating_dadi_SNP_input_from_structure#what-if-the-encoding-for-the-snps-is-0-3-for-a-t-and--9-for-missing for your purposes...

#How to run the script and what it spits out
