#################
##The user interface for our method
###################################
from os.path import isfile;
import pysnptools;
import sys;
from DP_util import *;


##
##Return the mret top scoring SNPs in an epsilon-DP way
##
def Top(MU,y,epsilon,mret,algor,savename):
    print "Choosing Top SNPS!"
    topSNPs=PickTop(y,MU,mret,epsilon,algor)

    print "The top "+str(mret)+" SNPs are estimated to be...";
    for s in topSNPs:
        print s;
    
    if len(savename)>0:
        fil=open(savename);
        fil.write("The results of looking for the top "+str(mret)+" SNPs with privacy parameter "+str(epsilon)+":\n");
        for s in topSNPs:
            fil.write(s+"\n");
        fil.close();
        print "SAVED!";
    print "\nGood Bye!"


##
##returns an estimate of the number of significant SNPs in an epsilon-DP way
##
def count(MU,y,epsilon,pval,savename):
    print "Counting SNPs!"
    
    num=estNum(MU,y,pval,epsilon);

    print "The estimated number of significant SNPs is: "+str(num);
    print "\nGood Bye!"

##
##returns an estimate of the Wald statistic on specified SNPs in an epsilon-DP way
##
def wald(MU,y,epsilon,snps,savename):
    print "Estimating Wald Statistics!"
    
    chi=WaldTest(y,MU,epsilon,snps);
    
    print "The scores are estimated to be..."
    for i in range(0,len(snps)):
        print snps[i]+" "+str(chi[i]);
    
    if len(savename)>0:
        fil=open(savename);
        fil.write("The results of calculating the Wald Statistic with privacy parameter "+str(epsilon)+":\n");
        for i in range(0,len(snps)):
            fil.write(snps[i]+" "+str(chi[i])+"\n");
        fil.close();
        print "SAVED!";
    print "\nGood Bye!"

##
##The user interface!
##
def Interface(args=[]):
	print "Hello, and welcome to PrivGWAS!"
	print "Checking dependencies, one moment..."
	try:
		import pysnptools;
	except ImportError, e:
		print "Sorry, you need pysnptools installed to run this program."
		return;
	try:
		import numpy;
	except ImportError, e:
		print "Sorry, you need numpy installed to run this program."
		return;
	try:
		import scipy;
	except ImportError, e:
		print "Sorry, you need scipy installed to run this program."
		return;
	try:
		import sklearn;
	except ImportError, e:
		print "Sorry, you need scikit learn installed to run this program."
		return;
	se2=.5;
	sg2=.5;
	epsilon=1.0;
	mret=3;
	pval=.05;
	snpsToEst=[];
	bedFil="";
	k=5;
	typ="Count";##Number count, return SNPs, estimate statistic
	meth="EigenStrat";
	algor="neighbor";
	useCov=False

	if len(args)>1:
		print "\nUnpacking arguments\n\n"
		args=args[1:];
		i=0;
		while i<len(args):
			if args[i][0]=="-":
				try:
					a=args[i];
					i=i+1;
					if i>len(args)-1:
						print "Sorry, incorrect argument!"
					elif a=="-q":
						if args[i]==1:
							useCov=True;
					elif a=="-se2":
						se2=float(args[i]);
					elif a=="-sg2":
						sg2=float(args[i]);
					elif a=="-k":
						k=int(args[i]);
					elif a=="-e":
						epsilon=float(args[i]);
					elif a=="-p":
						pval=float(args[i]);
						typ="Count"
					elif a=="-mret":
						mret=int(args[i]);
						typ="High"	
					elif a=="-s":
						snps=[];
						typ="Score"
						while i<len(args) and args[i][0]!="-":
							snps.append(args[i]);
							i=i+1;
						continue;
					elif a=="-bed":
						bedFil=args[i];
					elif a=="-meth":
						meth=args[i];
					elif a=="-a":
						algor=args[i];
					elif a=="-t":
						typ=args[i];
				except ValueError:
					print "Sorry, one of your arguments (the one corresponding to "+a+") is not of the correct type";
					return;
			i=i+1;
	if not isfile(bedFil+".bed"):
		print "Not a real BED file!"
		if ".bed" in bedFil:
			print "Note: do not include the .bed suffix in your filename"
			print "So if you had \"file.bed\", you should give the program the filename \"file\""
		return;
	if not epsilon>0:
		print "Epsilon needs to be at least zero"
		return;
	if not meth in ["EigenStrat","LMM"]:
		print "Not a valid method!" 
		return;
	if typ=="Count":
		if pval>1 or pval<0:
			print "pval must be between 0 and 1"
	if typ=="Top":
		if mret<1:
			print "mret needs to be >0"
	if typ=="Wald":
		if len(snps)==0:
			print "No SNPs provided!"
			return;

	
	print "As it stands: "
	print "Method: "+meth;
	print "BedFile: "+bedFil;
	print "epsilon: "+str(epsilon);
	print "Type: "+typ;
	if typ!="Count":
		print "Aglorithm: "+algor;
	if typ=="Top":
		print "mret: "+str(mret);
	if typ=="Wald":
		print "SNPs: "+str(snps);
	if typ=="Count":
		print "pvals: "+str(pval);
	print "\n\n\n";
	print "Load Data!"
	y=[];X=[];sFil=[];Q=[];
	[y,X,Q,sFil]=getData(bedFil,useCov);	
	print "Calculating MU matrix"
	MU=DP.getMU(y,X,Q=Q,se2=se2,sg2=sg2,k=k,meth=meth);	
	if typ=="Top":
		picks=DP.PickTopSNP(y,X,mret,epsilon=epsilon,k=k,se2=se2,sg2=sg2,MU=MU,meth=meth,algor=algor)
		if len(picks)==0:
			print "Bad argument!"
			return;
		print "The mret top scoring SNPs:"
		for i in picks:
			print sFil.sid[i];
			print sFil.pos[i]
	elif typ=="Wald":
		picks=sFil.sid_to_index(snps);#[snps.index(i) for i in snps];
		Scores=DP.estWald(y,X,epsilon,k=k,snps=picks,MU=MU,meth=meth,algor=algor);
		print "The estimated Wald scores are:";
		for i in range(0,len(picks)):
			print snps[i]+": "+str(Scores[i]);
	elif typ=="Count":
		print "Start counting algorithm!"
		nm=DP.estNum(MU,y,X,epsilon,k=k,se2=se2,sg2=sg2,MU=MU,pval=pval,meth=meth);
		print "There are about "+str(nm)+" SNPS with p-values less then "+str(pval); 
	
	else:
		print "Not a correct type of query!"
		return;
	print "Done!!";


if __name__=="__main__":
	args=sys.argv;
	if len(args)>1 and args[1][0]!="-":
		fil=open(args[1]);
		lines=fil.readlines();
		fil.close()
		args=[];
		args.append("UI.py")
		for l in lines:
			s=l.strip().split();
			args.extend(s);
	Interface(args);
