import os
import scipy.stats as ss
import sys
from add_exp import add_exp
import numpy as np
import anova
from random import shuffle

if len(sys.argv) != 6:
	print('python aseQtl.py [ase file] [date] [snps directory] [n] [ttest/utest]')
	sys.exit()

date = sys.argv[2]
cutoff = 1.0
#open ase file and read in
aseDict = add_exp(sys.argv[1])
if sys.argv[4] == "0":
	out = open("../results/"+sys.argv[2]+"/eqtl.ase.all","w")
else:
	out = open("../results/"+sys.argv[2]+"/permute/eqtl.ase.permute."+sys.argv[4],"w")
	cutoffFile = open('../results/'+sys.argv[2]+"/cutoff")
	cutoff = float(cutoffFile.readline().rstrip())

out.write('scaf	pac	locus	N af1	maf	ase_hom	ase_het	'+sys.argv[5][0]+'	p')

for file in os.listdir(sys.argv[3]):
	if file[-4:] != "snps":
		continue
	aseList = []
	snpFile = open(sys.argv[3]+file,'r')
	scaf = file.split('.')[1]
	pac = file.split('.')[0]

	indEnt = snpFile.readline().split()
	indList = indEnt[1:]

	
	if pac not in aseDict.keys(): #DO WE HAVE DATA FOR THIS?
		#print(pac)
		continue	

	sampleSize = 0
	aseLevels = aseDict[pac]

	for i in indList:
		if i not in aseLevels.keys(): #is there data for this individual?
			aseLevels[i] = "NA" #give this ind an "NA" for data
			continue
	if sys.argv[4] == "0":
		aseList = [aseLevels[i] for i in indList] 
	
	else: #permute the ase levels
	        myKeys = aseLevels.keys()
                myValues = aseLevels.values()
                shuffle(myValues)
                aseDictionary = dict(zip(myKeys,myValues))
		aseList = [aseDictionary[i] for i in indList] 

	for line in snpFile:
		snpEnt = line.split()
		snpName = snpEnt[0]
		
		genoDict = {'hom1':0, 'het':0,'hom2':0}
		myAseDict = {'hom':[],'het':[]}
		sampleSize = 0
		for i, genot in enumerate(snpEnt[1:]): 
			#add ase levels to genoDict
			if genot == "NA" or aseList[i] == "NA":
				continue
			else:
				sampleSize += 1
				genoDict[genot] += 1
				shortGen = genot[0:3]
				myAseDict[shortGen].append(float(aseList[i]))
		#do we have enough data?
		if len(myAseDict['het']) < 10 or len(myAseDict['hom']) < 10:
			print("not enough data "+snpName)
			continue
		
		#calculate allele frequency
		af1 = float(genoDict['hom1']+.5*genoDict['het'])/sampleSize

		#calculate minor allele freq
		if af1 > 0.5:
			maf = 1 - af1
		else:
			maf = af1

		#calc ase values for each category
		homMean = float(sum(myAseDict['hom']))/len(myAseDict['hom'])
		hetMean = float(sum(myAseDict['het']))/len(myAseDict['het'])

		#do the t test
		if sys.argv[5] == "ttest":
			myTtest = ss.ttest_ind(myAseDict['hom'],myAseDict['het'])
			tval, pval = myTtest[0],myTtest[1]
		elif sys.argv[5] == "utest":
			myUtest = ss.mannwhitneyu(myAseDict['hom'],myAseDict['het'])
			tval, pval = myUtest[0],myUtest[1]

		if sys.argv[4] == "0":
			out.write("\n"+"	".join([scaf,pac,snpName,str(sampleSize),str(af1),str(maf),str(homMean),str(hetMean),str(tval),str(pval)]))
		elif pval < cutoff:
			out.write("\n"+"	".join([scaf,pac,snpName,str(sampleSize),str(af1),str(maf),str(homMean),str(hetMean),str(tval),str(pval)]))


	snpFile.close()
