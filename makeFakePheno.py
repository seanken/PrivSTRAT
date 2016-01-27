#############
##Given genotype data
##creates fake phenotype data
##by picking one SNP and having the individual have P p0
##of having the disease if have 0 copie of minor allele,
##c*p0 if have 1, and c^2p0 if have 2
##where c is user defined, p0 chosen so about per percent of people have
##the disease
###################

import math;
import scipy as sp;
import numpy as np;
import random as rand;
from pysnptools.snpreader import Bed;

##
##Generates fakePheno based on data
##
def genPheno(filename="../thinFam",per=.5,savename="fakePheno.txt",c=2.0,num=5):
	sFil=Bed(filename);
	D=sFil.read().val;
	m=len(D[0]);
	n=len(D);
	print m;
	print n;
	I=[rand.randint(0,m-1) for i in range(0,num)];
	SNP=[[D[j][i] for j in range(0,n)] for i in I]
	#p0=n*peir/sum([c**i*len([j for j in SNP if j==float(i)]) for i in range(0,3)])
	print len(I);
	print len(SNP);
	print len(SNP[0]);
	print n;
	print min([len(s) for s in SNP])
	print SNP;
	
	SNP=[[max(i,0.0) for i in s] for s in SNP]
	for i in range(0,num):
		for j in range(0,n):
			if not SNP[i][j] in [1.0,0.0,2.0]:
				SNP[i][j]=0.0;
	print [list(set(s)) for s in SNP]
	lst=[sum([SNP[j][i] for j in range(0,num)]) for i in range(0,n)]
	#print lst;
	print sum([c**(sum([SNP[j][i] for j in range(0,num)])) for i in range(0,n)])
	p0=n*per/sum([c**(sum([SNP[j][i] for j in range(0,num)])) for i in range(0,n)])
	print p0;
	y=[float(rand.uniform(0,1)<p0*c**sum([SNP[j][i] for j in range(0,num)])) for i in range(0,n)]
	if len(savename)==0:
		return y;
	fil=open(savename,"w");
	for i in y:
		fil.write(str(i)+"\n");
	fil.close();


if __name__=="__main__":
	for c in [1.1,1.5,2.0]:
		genPheno(savename="fakePheno_"+str(c)+".txt",c=c);
	

	


