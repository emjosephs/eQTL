# eQTL

These are scripts for running a local eQTL analysis. 

I am still in the process of cleaning up and writing READMEs.

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

## Running an eQTL analysis: eQTL_kw_hets.py

This script can either run a basic eQTL analysis OR permute your data randomly and then run an eQTL analysis.

### regular eQTL usage

> python eQTL_kw_hets.py [expression file] [date] ht [genotype file directory] 0

the expression file is just a flat file table with your genotypes as columns and genes as rows. For example:

>pac     161A    163J    16A     80G     81J     83Z

>20885617        0.0551181102362 0.105263157895  0.112676056338  0.0809248554913 0.0901639344262 

>20885618        14.7244094488   7.00877192982   8.04929577465   9.53179190751   7.03278688525   


instead of date, you can use the folder name within your results folder that you would like your output file to be written in

### permuted eQTL usage

> python eQTL_kw_hets.py [exp file] [date] permuteAll [genotype file directory] [n]

before you run, create a folder called permute in the results folder for these outputs to be written to

n is the number of the permutation that you are running 

this script will output results from all tests. If you would like to save space by only outputting results that are below a certain significance threshold, create a text file in the dated results folder called cutOff with your cutoff value in it and run:

> python eQTL_kw_hets.py [exp file] [date] permute [genotype file directory] [n]


### output

the output files will be tables where each line is a test between SNP and gene expression the fields are:

1. scaffold
2. gene name
3. SNP/locus that was tested
4. frequency of the reference allele
5. minor, or folded, allele frequency
6. U statistic
7. P value
8. mean expression level of reference homozygote
9. mean expression of heterozygote
10. mean expression of alternate homozygote
11. frequencies of each genotype
12. hardy-weinberg deviation p value



## quantifying allele-specific expression with aseValuesByGene.py

This script takes in vcfs for RNA and DNA and outputs a table of allele-specific expression values (ie allelic imbalances) for all genes and all individuals using the RNA vcf information at informative heterozygous sites (genotype is inferred from the DNA vcf)

### Usage
>python aseValuesByGene.py [DNA vcf] [RNA vcf] [out file] [min snps per gene] [list of median expression levels]

min snps per gene is an integer. If set to >1, then it only returns and ASE value if there are more than the specified number of informative heterozygous snps in the gene.

the list of median expression levels is used to normalize the output by sequencing depths. These are median expression levels for all genes from a certain individual. This file should have two columns -- the first column is individual name and the second is median expression.

### Output
The output is a table where each row is a gene and each column is an individual and the values are allele-specific expression. 


