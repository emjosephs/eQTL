#takes in an eQTL file and outputs allele freq and deviation from HWE
#0	 1	  2	   3	  4	  5	 6	  7	         8	 9	 10	 11	 12	 13	  		14
#scaf    pac     locus   freq    fold    effect  dom     d_over_a        f       p       hom1    het     hom2    Rsquared        f(hom1/het/hom2)
#scaf8   20911662        9414108 0.801075268817  0.198924731183  1852.76344538   -2564.06193821  -1.38391220132  2.88186576785   0.0612160359523 3705.52689076   -711.298492832  0       0.0601869981805 59,31,3


import sys
from scipy import stats
import numpy

# hwe takes observed genotypes (in list) and the freq of ref allele (p), outputs p value
def hwe(obs, freq):
	hwObs = numpy.array([float(x) for x in obs])
	sampleSize = sum(hwObs)
	p = float(freq)
	hwExp = numpy.array([x*sampleSize for x in [p*p, 2*p*(1.0-p),(1.0-p)*(1.0-p)]])
	chisq = stats.chisquare(hwObs, hwExp, 1)
	return(chisq[1])

def __main__(): 

	eqtlFile = open(sys.argv[1],'r')
	out = open(sys.argv[1]+".hwe",'w')
	out.write("maf	obs	exp	p\n")

	for line in eqtlFile:
		ent = line.split()
		if ent[0] == "scaf":
			continue
	
		p = float(ent[3])
		maf = float(ent[4])
		#hwFreqs = [maf*maf, 2.0*maf*(1.0-maf), (1.0-maf)*(1.0-maf)]
		hwObs = numpy.array([float(x) for x in ent[14].split(',')])
		sampleSize = sum(hwObs)
		hwExp = numpy.array([x*sampleSize for x in [p*p, 2*p*(1.0-p),(1.0-p)*(1.0-p)]])


		chisq = stats.chisquare(hwObs, hwExp,1)

		out.write("	".join([str(maf),','.join([str(x) for x in hwObs]),','.join([str(x) for x in hwExp]),str(chisq[1]),"\n"]))




