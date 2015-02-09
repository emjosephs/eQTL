# reads in the htseq combined file and makes a dictionary of it, which it then returns
#dictionary key is the gene, value is itself a dictionary with key as individual and value as the expression level.
import sys


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


if __name__ == "__main__":
	print(add_exp('../data/htseq_norm/all.basic.norm'))







