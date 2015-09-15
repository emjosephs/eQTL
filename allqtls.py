# inputs a snp file and performs an anova
import os
import scipy.stats
import sys
import numpy
import random
import operator

def main():
	#generate expression dictionary
	expDict = add_exp(_args.exp_file)
	aseDict = add_exp(_args.ase_file)

	#randomize it if needed
	if _args.permute:
		permDic = makePermDic(_args.permute_no)
		permDic['150X'] = "150X"	
		#out = open('/cap1/emily.josephs/ase/results/'+_args.date+'/permute/eqtl.ase.permute.'+_args.permute_no,'w')
		out = open('/cap2/emily.josephs/ase/results/'+_args.date+'/permute/'+_args.outfile+'.permute.'+_args.permute_no,'w')
	
	else:
		#out = open('/cap1/emily.josephs/ase/results/'+_args.date+'/eqtl.ase.out.all','w')
		#out = open('/cap2/emily.josephs/ase/results/'+_args.date+'/eqtl.ase.out.all','w')
		out = open('/cap2/emily.josephs/ase/results/'+_args.date+'/'+_args.outfile,'w')
	out.write('scaf pac     locus   N.eqtl	N.ase	af1    maf	ase_hom	ase_het    h.ase       p.ase	hom1_mean	het_mean	hom2_mean	h.eqtl	p.eqtl')
	if _args.category_sizes:
		out.write("	hom1	hom2	het	hom")
	out.write("\n")

	for file in os.listdir(_args.snp_dir):   #skip non snps files -- could edit to only look at certain scaffolds or genes
        	if file[-4:] != 'snps':
                	continue
        	snpFile = open(_args.snp_dir+file,'r')

        	#find the scaffold
        	scaf = file.split('.')[1]

	        # read in the header of individuals
        	indLine = snpFile.readline()
       		indEntry = indLine.split()
        	pac, indList = str(indEntry[0]), indEntry[1:]
        	try:
                	expLevels = expDict[pac]
        	except: #it's a gene that was removed from analysis because of low exp or variation
                #	print("no exp: "+pac)
                	continue
		try:
			aseLevels = aseDict[pac]
		except:
		#	print("no ase: "+pac)
			continue
	        
		expLevels['150X'], aseLevels['150X'] = "NA","NA"
		#read in ase levels
        	if _args.permute: #we're permuting
			expList = [expLevels[permDic[i]] for i in indList]
			aseList = [aseLevels[permDic[i]] for i in indList]
        	else: #nope not permuting
			expList = [expLevels[i] for i in indList]	
			aseList = [aseLevels[i] for i in indList]	

		#read through snp lines
        	for line in snpFile:

			#do the eqtl analysis
			genoExpDict, snpName = makeGenoDict(line, expList)
			categorySizesExp = getLengths(genoExpDict) #returns a dictionary of type -- size

			af, fold, smallestCatSize, largestCatSize, comHomName, eN = calcFreqs(categorySizesExp)
			#note, if a gene has expression, there's data for every ind, so we can use the expression info to calculate allele freqs here.
	
        		if smallestCatSize < 10: #only does analysis if there's enough data for eqtl
                        	continue
			
			if _args.subsample: #subsample down to even out frequencies
				if largestCatSize > 40:
					subsampleDict1040(genoExpDict, [comHomName, 'het'])
				else:
					continue

			ehom1_mean, ehet_mean, ehom2_mean = getMeans(genoExpDict, ['hom1','het','hom2'])	
        	        ehstat,epval = scipy.stats.mstats.mannwhitneyu(genoExpDict[comHomName], genoExpDict['het'])
			#do the ase analysis
			genoAseDict, snpName = makeGenoDict(line, aseList)
        	        ahom_mean, ahet_mean = getMeans(genoAseDict,['hom','het'])
			
			catAse = getLengths(genoAseDict)
			if min(catAse['hom'], catAse['het']) < 10:
				ahstat,apval = "NA","NA"
        	        else:
				if _args.subsample:
					if max(catAse['het'], catAse['hom']) > 39:
					#subsampleDict(genoAseDict, ['hom','het'],10)
						subsampleDict1040(genoAseDict,['hom','het'])
					else:
						continue
				ahstat,apval = scipy.stats.mstats.mannwhitneyu(genoAseDict['hom'], genoAseDict['het'])
        	        	ahom_mean, ahet_mean = getMeans(genoAseDict,['hom','het'])
			
			aseN = str(catAse['hom']+catAse['het'])

                	outList = [scaf,pac,snpName,eN,aseN,af,fold, ahom_mean, ahet_mean, ahstat, apval, ehom1_mean, ehet_mean, ehom2_mean, ehstat, epval] 
			if _args.category_sizes:
				outList.extend([str(categorySizesExp[x]) for x in ['hom1','hom2','het','hom']])

                	out.write("     ".join([str(x) for x in outList])+"\n")

def subsampleDict(dic, namelist, n):
	"""
	>>> mydic = {"test":range(0,100), "ignore":range(0,100)}
	>>> subsampleDict(mydic, ['test'],10)
	>>> len(mydic["test"])
	10
	>>> len(mydic["ignore"])
	100
	"""
	for thing in namelist:
		if len(dic[thing]) < n:
			print('oops')
			sys.exit()
		#newList = [random.choice(dic[thing]) for x in range(0,n)]
		newList = random.sample(dic[thing], n)
		dic[thing] = newList
	
def subsampleDict1040(dic, namelist):
	"""
	>>> mydic = {"het":range(0,50), "hom":range(0,15)}
	>>> subsampleDict1040(mydic, ['het','hom'])
	>>> len(mydic["het"])
	40
	>>> len(mydic["hom"])
	10

	"""
	#print(dic)
	catSizes = dict(zip(namelist,[len(dic[x]) for x in namelist]))
	largerCat = max(catSizes.iteritems(), key=operator.itemgetter(1))[0]
	smallerCat = min(catSizes.iteritems(), key=operator.itemgetter(1))[0]
	newlargerCat = random.sample(dic[largerCat],40)
	dic[largerCat] = newlargerCat
	newsmallerCat = random.sample(dic[smallerCat],10)
	dic[smallerCat] = newsmallerCat
	


def makeGenoDict(line, expList):
	"""
	>>> line = "2736551 hom1    hom1    hom1"
	>>> expList = (3.0, 3.0, 3.0)
	>>> makeGenoDict(line, expList)
	({'het': [], 'hom': [3.0, 3.0, 3.0], 'hom1': [3.0, 3.0, 3.0], 'hom2': []}, '2736551')
	"""
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
	genoDict['hom'] = genoDict['hom1']+genoDict['hom2'] #combine hom data
	return(genoDict, snpName)

def myMean(myList):
        '''
        >>> list=[1.0,2,3]
        >>> myMean(list)
        2.0

        >>> empty=[]
        >>> myMean(empty)
        'NA'
        '''
        if(len(myList)) == 0:
                return('NA')
        else:
                return(sum(myList)/len(myList))

def getMeans(genoDict, tags):
	'''
	>>> genoDict = {'hom':[3.0, 3.0, 3.0], 'het':[3,0, 2.0, 1.0]}
	>>> getMeans(genoDict, ['hom','het'])
	[3.0, 1.5]
	'''	
	
	mymeans = [myMean(genoDict[x]) for x in tags]
	return(mymeans)


def getLengths(genoDict):
	'''
	>>> genoDict = {'hom':[3.0, 5.0, 6.0, 7.0],'hom1':[3.0], 'hom2':[5.0, 6.0, 7.0], 'het':[3.0, 2.0, 1.0, 1.0, 1.0]}
	>>> getLengths(genoDict)['hom'] 
	4
	'''
	names = ['hom1','het','hom2','hom']
	myLengths = [len(genoDict[x]) for x in names]
	return(dict(zip(names, myLengths)))

def calcFreqs(myDic):
        '''
	>>> calcFreqs({"hom1":25, "het":50, "hom2":25})
	(0.5, 0.5, 25, 50, 'hom2', 75)
	>>> calcFreqs({"hom1":0, "het":0, "hom2":0,"hom":0})
	(0, 0, 0, 0, 'NA', 0)
	>>> calcFreqs({"hom1":10, "het":40, "hom2":50,"hom":60})
	(0.3, 0.3, 40, 50, 'hom2', 90)
	>>> calcFreqs({"hom1":50, "het":40, "hom2":10,"hom":60})
	(0.7, 0.30000000000000004, 40, 50, 'hom1', 90)
	'''
	hom1, het, hom2 = [myDic[x] for x in ['hom1','het','hom2']]
	if hom1 + het + hom2 == 0: #there's no data here
		return(0,0,0,0, "NA",0)
	 
        af = float(2*hom1 + het)/(2*(hom1 + hom2 + het))
        if af <= 0.5:
        	fold = af #hom2 is common
		comHom, comHomName = hom2,"hom2"
        else:
                fold = 1 - af
		comHom, comHomName = hom1, "hom1"
	smallestCat = min([comHom, het])
	largestCat = max([comHom, het])
	eN = het+comHom
	return(af, fold, smallestCat, largestCat, comHomName, eN)

def makePermDic(n):
	'''
	>>> n='1'
	>>> makePermDic(n)['161A']
	'157x'
	>>> makePermDic('0')['161A']
	'161A'
	>>> makePermDic('2000')['161A']
	oops no match
	0
	'''
	permList = open('/cap1/emily.josephs/eQTL/bin/lists/paired_permuted_list','r')
        for line in permList:
        	sline = line.split()
                permNo = sline[0]
                if permNo == "0":
			realIndList = sline[1:]
		if permNo == str(n):
                	permIndList = sline[1:]
			return(dict(zip(realIndList, permIndList))) #have a dictionary of how to map permuted data
	print('oops no match')
	return(dict(zip(realIndList, [0]*len(realIndList))))

def add_exp(file):
	htDict = {}
	htFile = open(file, 'r')
        for line in htFile:
        	#get individual names
                if line[0:3] == "pac":
                        indList = line.split()
                else:
                        #get gene name
                        myEntry = line.split()
                        geneName = str(myEntry[0])
                        #get individual data
                        indDict = {}
                        for i, v in enumerate(myEntry):
                                if i == 0:
                                        continue
                                else:
                                        indName = indList[i]
                                        indDict[indName] = v
                        htDict[geneName] = indDict
        return(htDict)


import argparse
_args =  None
def parseArgs():
	parser = argparse.ArgumentParser(description="eqtl analysis")
	parser.add_argument("-e","--exp_file", type=str, help="expression table")
	parser.add_argument("-a","--ase_file", type=str, help="ase table")
	parser.add_argument("-d", "--date", type=str, help="date")
	parser.add_argument("-s", "--snp_dir", type=str, help="snp directory")
	parser.add_argument("-o", "--outfile", type=str, help="outfile name, default is eqtl.ase.out.all", default="eqtl.ase.out.all")
	parser.add_argument("-n", "--permute_no", type=str, help="permute no")
	parser.add_argument("-t", "--test", action="store_true", help="doctest?")
	parser.add_argument("-p", "--permute", action="store_true", help="permute the data?")
	parser.add_argument("-x", "--subsample", action="store_true", help="subsample down to 40/10 in each category")	
	parser.add_argument("-c", "--category_sizes", action="store_true", help="prints out the raw category sizes")	

	global _args
	_args = parser.parse_args()
	sys.stderr.write(str(_args)+"\n")


if __name__ == "__main__":
	parseArgs()

	if _args.test:
		import doctest
		doctest.testmod()
		sys.exit()
	main()


		
