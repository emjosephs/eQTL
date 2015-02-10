# eQTL

These are scripts for running a local eQTL analysis. 

I am still in the process of cleaning up and writing READMEs.

The main idea is that you want to keep these scripts in a bin folder and create a "data" folder and a "results" folder adjacent to your scripts folder. Inside the results folder, mmake specific folders for each run -- I do this using dates.

Here are the steps for analysis:
## Creating genotype files for each gene: make_snps_files.py

### Usage
> python make_snps_files.py [gff file] [vcf file] [scaf_number] [distance] [out directory] [all/noncoding] [annotation summary file]

You need [PyVCF](https://github.com/jamescasbon/PyVCF) to run this script.

Distance refers to the distance from the TSS or TES that you want to include in the file.

With the all option, you get all SNPs included in the genotype file and with the noncoding option you only get noncoding SNPs If you use the noncoding option, you need to add an annotation summary file. If not, you can use a "." here. [You can read more about the annotation format here](http://www.genomicconflict.com/wiki/index.php?title=Roberts_Annotations)

### Output
After running this  you will end up with a bunch of files in your out directory that are named [gene name].[scaffold].snps 

For example: 20889218.scaf1.snps

The first line of this file will have the gene name and all the individuals. The subsequent lines will have the SNP coordinate and then genotypes at this coordinate.

They should look like this:

>20889218	16A	11G	138Q	203A	114I	24F
>
>1577995	het	het	het	hom1	hom1	hom1
>
>1578031	hom2	hom2	het	het	het	hom2
>
>1578033	hom1	hom1	hom1	hom1	hom1	hom1

note that hom1 is the reference homozygote

## Running an eQTL analysis

