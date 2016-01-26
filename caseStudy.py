from testCase import loadData;
from MU_STRAT import *;
from DP_util import WaldTest as wt;
from DP_util import CI;
from DP_util import PickTop as pt;


def TestCase(cutoff=23,num=2,filename="../GWAS/NoFam_0.05",k=10,sz=450):
	[[y1,sFil1],[y2,sFil2]]=loadData(filename=filename,useRand=False,sz=sz);
	MU1=MU_STRAT(sFil1,k)
	MU2=MU_STRAT(sFil2,k)

	n=len(y2);
	n2=len(y1);
	print "Exact 1!"
	MU1.calcMU(k,exact=True);
	print "Exact 2!"
	MU2.calcMU(k,exact=True);

	print "Get top in set 1"
	top1=pt(y1,MU1,20,-1,algor="noise");
	print "Get Scores"
	scores=wt(y1,MU1,-1,snps=top1);

	print "Get "+str(num)+" SNPs directly below "+str(cutoff);
	i=min([i for i in range(0,len(top1)) if scores[i]<cutoff/float(n2-k-1)])
	print top1[i];
	print float((n2-k-1)*scores[i])
	#print i;
	#print scores[:10];
	print "";
	print top1[i+1];
	print float((n2-k-1)*scores[i+1])

	CIlst=CI(y2,MU2,.95,2.0,top1[i:i+2]);
	Clst=[[(n-k-1)*i for i in s] for s in CIlst]
	print Clst;




if __name__=="__main__":
	TestCase(cutoff=24,num=2,filename="../GWAS/NoFam"k=10,sz=450);




