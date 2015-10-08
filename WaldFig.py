####################
##Generates figure similar to those in the paper
#####################
import sys;
import loadFile as ld;
from MU_STRAT import *;
from DP_util import WaldTest as wt;



##
##For making pics!
##
def plotWald(eps,filename,savename="",k=5):
	if len(savename)==0:
		savename="OutputDir/res_wald_"+str(eps[0])+"_"+str(k)+".txt"
	print "load Data"
	[y,BED]=ld.getData(filename);
	print "calc MU!"
	MU=MU_STRAT(BED,k);
	print "get True!";
	tru=wt(y,MU,-1,snps=[],forFigs=False);
	fil=open(savename,"w");
	n=len(y);
	for i in range(0,10):
		e=eps[i];
		print e;
		res=wt(y,MU,e,snps=[],forFigs=True);
		err=sorted([float(n-k-1)*abs(res[i]-tru[i]) for i in range(0,len(tru))]);
		m=len(err);
		med=err[int(.5*m)]
		up=err[int(.75*m)];
		down=err[int(.25*m)];
		print med;
		fil.write(str(e)+" "+str(down)+" "+str(med)+" "+str(up)+"\n");
	fil.close();



if __name__=="__main__":
	filename="../../GWAS/cleaned";
	eps=[.5*i for i in range(1,11)];
	savename="OutputDir/cleaned_Wald";

	plotWald(eps,filename,savename,k=5);
	
	filename="../FakeData/pop"
	eps=[.5*i for i in range(1,11)];
	savename="OutputDir/pop_Wald";

	plotWald(eps,filename,savename,k=5);
