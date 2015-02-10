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

After running this  you will end up with a bunch of files in [out directory] that are named [gene name].[scaffold].snps 

The first line of this file will have the gene name and all the individuals. The subsequent lines will have the SNP coordinate and then genotypes at this coordinate.

They should look like this:

>20889218	16A	11G	138Q	203A	114I	24F	8U	129Y	50S	155T	121B	123D	144Q	58Y	106A	63D	124Y	175D	147Y	195L	194F197L	190U	193N	192H	115L	89L	111C	110M	113A	112H	83Z	86I	54T	143Y	202I	80G	109E	163J	103K	100O	101S	107B	139G102K	34D	150X	44W	49Z	61G	152P	67R	177K	64E	154P	126Z	174I	183E	181Y	187K	189B	85I	23Z	91C	93M	92H	95O94U	96O	10M	65T	78M	151D	153G	81J	97C	157x	158J	170V	1N	14R	98S	55Y	99X	9C	7K	146J	200U	39N	140S117R	74M	72G	71W	137Q	79G	161A	43F	4S	125X
>1577995	het	het	het	hom1	hom1	hom1	hom1	hom1	hom1	hom1	NA	hom1	het	hom1	hom1	hom1	het	hom1	hom1	het	hom1	hom1het	het	hom1	hom1	hom1	hom1	het	het	hom1	het	hom1	hom1	hom1	hom1	hom1	hom1	het	hom1	hom1	hom1	hom1	hom1	hom1het	hom1	hom1	hom1	NA	hom1	hom1	hom1	hom1	hom1	hom1	hom1	het	hom1	hom1	het	hom1	hom1	hom1	hom1	hom1	hom1	hethom1	hom1	hom1	hom1	hom1	NA	het	hom1	hom1	hom1	hom1	hom1	hom1	het	hom1	het	hom1	hom1	hom1	hom1	het	hom1	hethet	het	hom1	hom1	het	hom1	hom2	hom1	hom1
>1578031	hom2	hom2	het	het	het	hom2	het	het	hom2	hom2	het	het	hom2	hom2	het	het	hom2	het	hom2	hom2	hom2	hom2hom2	NA	hom2	hom2	hom2	hom2	hom2	hom2	hom2	hom2	hom1	hom2	hom2	het	het	het	hom2	hom2	hom2	hom2	het	hom2	hethom2	hom2	hom2	hom2	NA	hom2	hom2	het	het	hom2	hom2	het	hom2	het	het	hom2	het	hom2	hom2	het	het	hom1	hethom2	het	hom2	hom2	het	het	hom2	hom2	hom2	hom2	hom2	hom2	hom2	hom2	hom2	hom2	hom2	hom2	hom2	hom2	hom2	het	hom2hom2	hom2	hom2	hom2	hom2	het	hom2	het	hom2
>1578033	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	NA	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1hom1	NA	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1hom1	hom1	hom1	hom1	NA	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1hom1	hom1	hom1	hom1	hom1	NA	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1	hom1



