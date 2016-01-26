##
##Demonstrates test case
##
import matplotlib.pyplot as plt;
from MU_STRAT import *;
from DP_util import WaldTest as wt;
from loadFile import getData;
from DP_util import PickTop as pt;
from DP_util import CI;
import random as rand;
##
##loads both dataset 1 and dataset 2
##
def loadData(filename="",sz=100,useRand=False):
	if len(filename)==0:
		filename="../../GWAS/cleaned";
	[y,sFil]=getData(filename);
	if useRand:
		rand.shuffle(y);

	I=[];
	I1=[i for i in range(0,len(y)) if y[i]==1];
	I2=[i for i in range(0,len(y)) if y[i]!=1];
	
	I.extend(I1[:sz]);
	I.extend(I2[:sz]);

	NI=[];
	NI.extend(I1[sz:])
	NI.extend(I2[sz:])

	y1=[y[i] for i in I]
	sFil1=sFil[I,:];


	y2=[y[i] for i in NI]
	sFil2=sFil[NI,:];

	return [[y1,sFil1],[y2,sFil2]];


##
##Finds highest scoring SNP in database 1, returns, than queries 
##database 2 in a DP manner
##
def RunGWAS(filename,k=5,eps=2.0):
	sz=450
	[[y1,sFil1],[y2,sFil2]]=loadData(filename=filename,useRand=False,sz=sz);

	MU1=MU_STRAT(sFil1,k);	
	MU2=MU_STRAT(sFil2,k);
	print "EXACT1!"
	MU1.calcMU(k,exact=True)
	print "Exact 2!"
	MU2.calcMU(k,exact=True);

	res=pt(y1,MU1,1,-1,algor="noise");

	resL=pt(y1,MU1,20,-1,algor="noise");
	print resL
	#print wt(y1,MU1,-1,snps=res);
	res2=pt(y2,MU2,1,-1,algor="noise");
	print res;
	resT=res;
	print res2;
	resT=wt(y2,MU2,eps,snps=resL[5:7]);
	res2=wt(y1,MU1,eps,snps=res);
	res3=wt(y1,MU1,-1,snps=[]);

	res4=wt(y2,MU2,-1,snps=[]);
	n=len(y2);

	print len(y1);
	print n;
	val=[resT[0]*(n-k-1),resT[1]*(n-k-1)]

	n2=len(y1);	
	print "In Val:"
	print val;
	print "In Orig"
	print res2[0]*(n2-k-1);

	n2=len(y1);	
	print np.median(res3)*(n2-k-1)

	n2=len(y1);	
	res3=sorted(res3,reverse=True)[:10];
	res3=[(n2-k-1)*i for i in res3];
	print res3

	fil=open("saveRes.txt","w")
	for r in res4:
		fil.write(str(r)+"\n")
	fil.close();
	res4=sorted(res4,reverse=True)[:20];
	res3=[(n-k-1)*i for i in res4];
	print res3
	CIlst=CI(y2,MU2,.95,2.0,resL[5:6]);
	CIlst2=CI(y2,MU2,.95,2.0,resL[6:7]);
	CIlst=[[(n-k-1)*i for i in s] for s in CIlst]
	CIlst22=[[(n-k-1)*i for i in s] for s in CIlst2]
	print CIlst;
	print CIlst2;

if __name__=="__main__":
	filename="../GWAS/NoFam_0.05"
	#filename="../../HapAndRA/HapAndRA2"
	RunGWAS(filename,k=5,eps=2.0);
