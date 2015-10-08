######################
##Loads Plink file
######################
import math;
from numpy.random import permutation as perm
import numpy as np;
import Useful as us;
import os;
from pysnptools.util import intersect_apply;
from pysnptools.snpreader import Bed;
from pysnptools.snpreader import Pheno;
from os.path import isfile;

##
##gets y, the phenotype, X, the normalized matrix, sFil,
##the snpreader, and Q the covariate matrix
##
def getData(filename="",mph=3,UseCov=False):
	sFil=Bed(filename);
	yFil=Pheno(filename+".fam");
	
	Q=[];
	if isfile(filename+".cov") and UseCov:
		QFil=Pheno(filename+".cov")
		[sFil,yFil,QFil]=intersect_apply([sFil,yFil,QFil])				
	if isfile(filename+".phen"):
		yFil=Pheno(filename+".phen");
		[sFil,yFil]=intersect_apply([sFil,yFil])				
	return [yFil,sFil];
	

##
##Calculates se2 and sg2 for given filename
##
def varRes(filename,direct):
	print "Estimating GRM"
	os.system("gcta64 --bfile "+filename+" --make-grm --autosome --out "+filename+" > "+direct+"/tmp.txt")
	print "Calculating Heritability"
	os.system("gcta64 --reml --grm "+filename+" --pheno "+filename+".phen --out "+filename+" > "+direct+"/tmp.txt");
	os.system("rm "+direct+"/tmp.txt")
	if isfile(filename+".hsq"):
		fil=open(filename+".hsq");
		lines=fil.readlines();
		fil.close();
		se2=0.0;
		sg2=0.0;
		for l in lines:
			s=l.split();
			if s[0]=="V(G)":
				sg2=float(s[1]);
			if s[0]=="V(e)":
				se2=float(s[1]);
		return [se2,sg2];
	return [-1.0,-1.0];

##
##Divide data!
##
def divideData(filename,direct,num=5,mph=3,delet=True):
	print "Estimating heritability using "+str(num)+" components"
	[yFil,sFil]=getData(filename,mph=mph);
	n=sFil.iid_count	
	reOrd=perm(n);
	yFil=yFil[reOrd,:];
	sFil=sFil[reOrd,:];

	div=[int(math.ceil( i*n/float(num) )) for i in range(0,num+1)];
		
	varEsts=[];

	for i in range(0,num):
		print "For component "+str(i);
		sFilTemp=sFil[div[i]:div[i+1],:];

		yFilTemp=yFil[div[i]:div[i+1],:];

		fileTemp=direct+"/tempFile_"+str(i);
		Bed.write(fileTemp,sFilTemp.read());
		Pheno.write(fileTemp+".phen",yFilTemp.read())
		
		varEsts.append(varRes(fileTemp,direct));
		
		

		if delet:
			os.system("rm "+direct+"/tempFile_"+str(i)+"*");
	
	return varEsts;


##
##estimates heritability
##
def estHerit(filename,num,epsilon,getVar=False,mph=3,direct="Temp"):
	varEst=divideData(filename,direct,num=num,mph=mph,delet=True);
	if getVar:
		return varEst[0];
	est=sum([v[1]/(sum(v)) for v in varEst])/float(num);
	if epsilon<0:
		return [sum(varEst[0]),est];
	print "Averaging and adding noise"	
	est=est+us.Lap(0,epsilon/float(num));
	if est<0:
		est=0.0;
	if est>1:
		est=1.0
	print "The Result is "+str(est);
	return est;

##
##estimates the variance of y
##
def estVarY(y,epsilon):
	vr=np.var(y);
	n=len(y);
	return vr+us.Lap(0,epsilon*3/float(n))



if __name__=="__main__":
	filename="../Nice_DP_Pop/FakeData/pop"
	h2=[];
	h1=estHerit(filename,num=1,epsilon=-1)
	for i in range(0,1):
		h2.append(estHerit(filename,num=5,epsilon=1.0))
	print h1;
	print h2;


	#divideData("../GWAS/cleaned","Temp",num=1,delet=True)
