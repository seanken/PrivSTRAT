#################
##The user interface for our method
###################################
from os.path import isfile;
import pysnptools;
import sys;
from loadFile import getData;
from MU_LMM import MU_LMM;
import PrivGWAS;

##
##The user interface!
##
def Interface(args=[]):
	print "Hello, and welcome to PrivLMM!"
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
	epsilon=1.0;
	mret=3;
	pval=.05;
	snpsToEst=[];
	bedFil="";
	num=0;
	se2=.5;
	sg2=.5;
	typ="Top";##Number count, return SNPs, estimate statistic
	algor="noise";
    	savename="";


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
					if a=="-se2":
						se2=float(args[i]);
					if a=="-sg2":
						sg2=float(args[i])
					if a=="-num":
						num=int(args[i])
					if a=="-e":
						epsilon=float(args[i]);
					if a=="-p":
						pval=float(args[i]);
						typ="Count"
					if a=="-mret":
						mret=int(args[i]);
						typ="Top"	
					if a=="-s":
						snps=[];
						typ="Wald"
						while i<len(args) and args[i][0]!="-":
							snps.append(args[i]);
							i=i+1;
						continue;
					if a=="-bed":
						bedFil=args[i];
        		            	if a=="-save":
                        			savename=args[i];
					if a=="-a":
						algor=args[i];
					if a=="-t":
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
	if typ=="Herit":
		if num<1:
			print "num must be > 0";
			return;
		se2=-1.0;
		sg2=-1.0
	
	print "As it stands: "
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
    	[y,BED]=getData(bedFil);

	print "Calculating MU matrix"
	if num<1:
		MU=MU_LMM(BED,(se2,sg2));
	else:	
		MU=MU_LMM(BED,[num,epsilon]);
	if typ=="Top":
        	PrivGWAS.Top(MU,y,epsilon,mret,algor,savename);
   	elif typ=="Count":
        	PrivGWAS.count(MU,y,epsilon,pval,savename);
	elif typ=="Wald":
        	PrivGWAS.wald(MU,y,epsilon,snps,savename);
	elif typ=="Herit":
		print "The estimated heritability:"
		print "Sigma_e^2 is "+str(MU.se2);
		print "Sigma_g^2 is "+str(MU.sg2);
		print "Good bye!"
	else:
		print "Specified task "+typ+" not recognized."
		print "Good bye!"



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
