# eQTL

These are scripts for running a local eQTL analysis. 

I am still in the process of cleaning up and writing READMEs.

The main idea is that you want to keep these scripts in a bin folder and create a "data" folder and a "results" folder adjacent to your scripts folder. Inside the results folder, mmake specific folders for each run -- I do this using dates.

Here are the steps for analysis:
## Creating genotype files for each gene: make_snps_files.py

usage: python make_snps_files.py [gff file] [vcf file] [scaf_number] [distance] [out directory] [all/noncoding] [annotation summary file]

Distance refers to the distance from the TSS or TES that you want to include in the file.

With the all option, you get all SNPs included in the genotype file and with the noncoding option you only get noncoding SNPs

If you use the noncoding option, you need to add an annotation summary file. If not, you can use a "." here.

[You can read more about the annotation format here](http://www.genomicconflict.com/wiki/index.php?title=Roberts_Annotations)




