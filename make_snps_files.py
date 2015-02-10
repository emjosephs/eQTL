# this is actually the master script -- bad name oops


import vcf
import sys
import gff_parse
from name_switch import makeDict
import pickle
import vcf_parse
from add_exp import add_exp


def __main__():

	#check arguments
	if len(sys.argv) != 8:
		print(len(sys.argv))
		print('python make_snps_files.py [gff file] [vcf file] [scaf_number] [distance] [out directory] [all/noncoding]')
		sys.exit()

	#nameDict = makeDict('ind_list_4-2-2014')
	nameDict= makeDict('ind_lists/ind_list_4-21-2014_nofield')
	nameList = []
	for ind in nameDict.values():
		nameList.append(ind)
	print(nameList)

	# read in the gff
	parseOut = gff_parse.gff_parse(sys.argv[1],int(sys.argv[4]))
	gffList = parseOut['gffList']
	geneDict = parseOut['geneDict']

	#read in the annotation file
	if sys.argv[6] == 'noncoding':
		annot = open(sys.argv[7],'r')
		annotDic =  makeAnnotDic(annot)

	# read in the vcf
	vcf_reader = vcf.VCFReader( open(sys.argv[2],'r') )
	for record in vcf_reader:
		
		entry = vcf_parse.vcfParse(record, nameDict)
		# see if it's in a gene!
		if len(gffList[entry['pos']]) == 0:
			continue
		# is it in a coding site and are we making noncoding?
		elif sys.argv[6] == 'noncoding':
			if annotDic[entry['pos']] == 1: #make it an integer
				continue
		#add the snp to each gene
		for gene in gffList[entry['pos']]:
			geneDict[gene].append(entry)
	#print(geneDict['20911598'][0])
	
	# write out.
	for gene in geneDict:
		out = open(sys.argv[5]+gene+".scaf"+str(sys.argv[3])+".snps","w")
		
		#write out first line of individual names
		out.write(gene)
		for ind in nameList:
			out.write("	"+str(ind))

		#write out snp genotypes
		for snp in geneDict[gene]:
			out.write("\n"+str(snp['pos']))
			#print(geneDict[gene])
			for ind in nameList:
				out.write("	"+str(snp["genotypes"][ind]))


def makeAnnotDic(annot): #makes a dictionary of all sites, 0 if noncoding, 1 if in an exon
	exonIDs = []
        annotDic = [0]*20000000
	#exonIDs = [3,2,4]
        for line in annot:
		ent = line.split()
                if line[0:5] == "#TYPE" and ent[1] in ['0fold','exon','4fold','stop']:
                        exonIDs.append(ent[2].rstrip())
                elif line[0] == "#":
                        pass
                else:
                        scaf,pos,sts = int(ent[0][-1]),int(ent[1]), ent[7]
			if any(x in exonIDs for x in sts.split(',')): #this is an exon site
        			annotDic[pos] = 1
	return(annotDic)



if __name__ == "__main__":
	__main__()

