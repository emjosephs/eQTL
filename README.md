# eQTL

These are scripts for running a local eQTL and aseQTL analysis. 

The main idea is that you want to keep these scripts in a bin folder and create a "data" folder and a "results" folder adjacent to your scripts folder. Inside the results folder, mmake specific folders for each run -- I do this using dates.

To use these you will need to install SciPy and NumPy

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


## quantifying allele-specific expression with aseValuesByGene.py

This script takes in vcfs for RNA and DNA and outputs a table of allele-specific expression values (ie allelic imbalances) for all genes and all individuals using the RNA vcf information at informative heterozygous sites (genotype is inferred from the DNA vcf)

### Usage
>python aseValuesByGene.py [DNA vcf] [RNA vcf] [out file] [min snps per gene] [list of median expression levels]

min snps per gene is an integer. If set to >1, then it only returns and ASE value if there are more than the specified number of informative heterozygous snps in the gene.

the list of median expression levels is used to normalize the output by sequencing depths. These are median expression levels for all genes from a certain individual. This file should have two columns -- the first column is individual name and the second is median expression.

### Output
The output is a table where each row is a gene and each column is an individual and the values are allele-specific expression. 


## Running eQTL and aseQTL analyses
This combined script will run both eQTL and aseQTL analysis on a set of snp files (generated above). 

Basic usage is:
> python allqtl.py -e [EXP_FILE] -a [ASE_FILE] -s [SNP DIRECTORY] -s [OUTPUT FILE]

If you would like to permute the data randomly

> python allqtl.py -e [EXP_FILE] -a [ASE_FILE] -s [SNP DIRECTORY] -s [OUTPUT FILE] -p -n [PERMUTATION NUMBER]

To subsample so that analyses are run on samples of 50 individuals (10 and 40 in each category).

> python allqtl.py -e [EXP_FILE] -a [ASE_FILE] -s [SNP DIRECTORY] -s [OUTPUT FILE] -x

To see all options:
> python allqtls.py -h

The output is a table with information for each snp. The columns are
scaf: scaffold name 
pac: gene name     
locus: SNP location   
N.eqtl: The number of individuals that could be tested for eQTLs	
N.ase: The number of individuals that could be tested for aseQTL	
af1: Frequency of the reference allele    
maf: Frequency of the minor allele	
ase_hom: Mean ASE of homozygotes for the SNP	
ase_het: Mean ase for heterozygotes for the snp    
h.ase: Mann-Witney statistic for aseQTL test       
p.ase: P value for aseQTL test	
hom1_mean: Mean expression of individuals homozygous for reference allele	
het_mean: Mean expression of individuals heterozygous at SNP	
hom2_mean: Mean expression of individuals homozygous for alternate allele	
h.eqtl: Mann-Whitney statistic for eQTL test	
p.eqtl: P value for eQTL test








