# inputs a snp file and performs an anova
import os
import scipy.stats
from add_exp import add_exp
import sys
from random import shuffle
import numpy
from calc_HW_deviations import hwe

def myMean(myList):
	if(len(myList)) == 0:
		return('NA')
	else:
		return(sum(myList)/len(myList))

if len(sys.argv) < 5:
	print('python eQTL_kw_hets.py [exp file] [date] [ht/permute/permuteAll] [snps directory] [permute no]')
	sys.exit()

date = sys.argv[2]

if sys.argv[3] == "ht":
	out = open('../results/'+date+'/eqtl.kw.out.all','w')
elif sys.argv[3] in ["permute","permuteAll"]:
	out = open('../results/'+date+'/permute/eqtl.kw.permute.'+sys.argv[5],'w')
	cutoffFile = open('../results/'+date+"/cutoff",'r')
	cutoff = float(cutoffFile.readline().rstrip())
err = open('../results/'+date+'/eqtl.anova.err','w')
out.write('scaf	pac	locus	freq	fold	h	p	hom1	het	hom2	f(hom1/het/hom2)	hwe_dev\n')
#generate expression dictionary
expDict = add_exp(sys.argv[1])

for file in os.listdir(sys.argv[4]):
	#skip non snps files -- could edit to only look at certain scaffolds or genes
	if file[-4:] != 'snps':
		continue
	htList = []
	snpFile = open(sys.argv[4]+file,'r')

	#find the scaffold
	scaf = file.split('.')[1]

	# read in the header of individuals
	indLine = snpFile.readline()
	indEntry = indLine.split()
	pac, indList = str(indEntry[0]), indEntry[1:]
	try:
		expLevels = expDict[pac]
	except: #it's a gene that was removed from analysis because of low exp or variation
		err.write("no expression: "+pac+"\n")
		continue	

	#read in expression levels from ht_dict
	if sys.argv[3] =="ht" or int(sys.argv[5]) == 0: #if 0, prints out permute thing for the real data!!
		pass
	elif sys.argv[3] in ['permute', 'permuteAll'] and int(sys.argv[5]) > 0:
		myKeys = expLevels.keys()
		myValues = expLevels.values()
		shuffle(myValues)
		expLevels = dict(zip(myKeys,myValues))	
	sampleSize = 0

	for i in indList:
		try:
			expLevels[i]
			sampleSize = sampleSize + 1
		except:
			expLevels[i] = "NA"
			err.write(i+"\n")
	expList = [expLevels[i] for i in indList]
	#err.write("sample size is "+str(sampleSize)+"\n")

	# read through snp lines
	for line in snpFile:
		snpEntry = line.split()
		snpName = str(snpEntry[0])
		
		#make bins for expression lebel
		genoDict = {'hom1':[],'het':[],'hom2':[]}

		#read through genotypes and end expression level to write bins
		for pos, ind in enumerate(snpEntry[1:]):
			if ind == 'NA' or expList[pos] == "NA":
				continue
			else:
				genoDict[ind].append(float(expList[pos]))

		#calculate the allele frequency of allele 1 
		if len(genoDict['het']+genoDict['hom2']+genoDict['hom1']) > 0:
			af = float(2*len(genoDict['hom1']) + len(genoDict['het']))/(2*(len(genoDict['hom1'])+len(genoDict['hom2'])+len(genoDict['het'])))
		else:
			err.write("no data: "+line)

		#calculate fold!
		if af <= 0.5:
			fold = af
			#hom2 is common
			comHom = genoDict['hom2']
		else:
			fold = 1 - af
			comHom = genoDict['hom1']

		#do we have enoughdata?
		if len(comHom) < 10 or len(genoDict['het']) < 10:
			continue

		#just publish raw means
		hom1_mean = myMean(genoDict['hom1'])	
		hom2_mean = myMean(genoDict['hom2'])	
		het_mean = myMean(genoDict['het'])		
		
		hwObs = [len(genoDict['hom1']),len(genoDict['het']),len(genoDict['hom2'])]
		dev = hwe(hwObs, af)
		gf = str(len(genoDict['hom1']))+","+str(len(genoDict['het']))+","+str(len(genoDict['hom2']))
		if dev < 0.05:
			err.write("hwe dev:	"+str(af)+"	"+gf)
			continue

		#do the kw test!
		hstat,pval = scipy.stats.mstats.mannwhitneyu(comHom, genoDict['het'])		
#		try:
#			hstat, pval = scipy.stats.mstats.kruskalwallis(comHom,genoDict['het'])
#		except ValueError:
#			err.write("ValueError	"+pac+"\n")
#			continue


		if sys.argv[3] in ["ht","permuteAll"]:
			out.write("	".join([scaf,pac,snpName,str(af),str(fold),str(hstat),str(pval),str(hom1_mean),str(het_mean),str(hom2_mean),gf,str(dev)])+"\n")
		
		elif sys.argv[3] == "permute" and pval != "nan" and pval < cutoff:
			out.write("	".join([scaf,pac,snpName,str(af),str(fold),str(hstat),str(pval),str(hom1_mean),str(het_mean),str(hom2_mean),gf,str(dev)])+"\n")

