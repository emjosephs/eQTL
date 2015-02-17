# goal = make a table of ase values for each gene.
from name_switch import makeDict
import sys
import scipy.stats
import vcf

if len(sys.argv) != 6:
	print('python aseValuesByGene.py [DNA vcf] [RNA vcf] [out file] [min snps per gene] [list of median expression levels]')
	sys.exit()

nameDict= makeDict('/cap1/emily.josephs/eQTL/bin/ind_lists/ind_list_5-13-2014_nofield')
dnaDict = {} #keys are (scaffold,pos) and value is a list of all the hets.
aseDict = {} #keys are (gene) value is dict where keys = ind, vals = list of diffs in het snps
med = {} #keys are ebjInd, values are the median.
indList = nameDict.values()
#make med
medList = open(sys.argv[5],'r')
for line in medList:
        ent = line.split()
        med[ent[0]] = float(ent[1])

#read in DNA genotypes and add to dnaDict
vcf_reader = vcf.Reader( open(sys.argv[1],'r') )
for record in vcf_reader:
        hets = [x.sample for x in record.get_hets()]
        dnaDict[(record.CHROM,record.POS)] = hets

#read through RNA genotypes and add to aseDict
rna_reader = vcf.VCFReader( open(sys.argv[2],'r'))
for record in rna_reader:
        
        if record.CHROM not in aseDict.keys():
                aseDict[record.CHROM] = {} #add gene to aseDict if it's not already there
		for i in indList:
			aseDict[record.CHROM][i] = []	

	for ind in dnaDict[(record.CHROM,record.POS)]: #read through inds with het genotypes
		#is this ind in our 99 samples?
		if ind not in nameDict.keys():
			continue
		ebjInd = nameDict[ind] #switch YW id to EBJ id
		
		if ebjInd not in aseDict[record.CHROM].keys():
			aseDict[record.CHROM][ebjInd] = [] #add gene to aseDict keys if not there already
		try:
			ad = record.genotype(ebjInd).data.AD
		except:
			print(ebjInd)
			print(record)
			print(record.samples)
			print('')  #want to make sure I'm not missing cases of complete ASE
			continue
		if ad == None: #missing data
			continue

		aseDif = float(abs(ad[1] - ad[0]))/med[ebjInd]
		aseDict[record.CHROM][ebjInd].append(aseDif)

# print out --  table where rows are genes, columns are genotypes
out = open(sys.argv[3],'a')
#out.write("pac	"+"	".join(indList))

for gene in aseDict:
	out.write(gene)
	for ind in indList:
		if len(aseDict[gene][ind]) >= int(sys.argv[4]): #do we have data?
			aseVal = sum(aseDict[gene][ind])/len(aseDict[gene][ind]) #calc aseMean for this ind
			out.write('	'+str(aseVal))
		else:
			out.write('	NA')
	out.write("\n")
