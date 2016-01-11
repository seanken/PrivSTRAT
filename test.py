##########
##implements some basic tests
#####################

import random as rand;
from MU_Mat import *
from DP_util import *;
from loadFile import *;
from MU_STRAT import *;
from MU_LMM import *;
from PrivSTRAT import *;
from pysnptools.snpreader import Bed; 

##
##Makes sure loads file correctly
##
def testloadfile():	
	print "Test load file!"
	passTests=True;
	testfile1="TestCases/pop";
	try:
		[y,sFil]=getData(testfile1);
	except:
		print "Exception thrown loading genotype data!"
		print "NOT PASS!!!\n\n\n"
		return;
	print "Loaded data!";
	if not isinstance(sFil,pysnptools.snpreader.Bed):
		print "returned value is not snpreader!"
		passTests=False;
	else:
		print "returns snpreader.Bed object!"
	if sum(y)==5000:
		print "Correct number of cases!"
	else:
		passTests=False
		print "Error loading phenotype!"
	if sFil.iid_count==10000 and sFil.sid_count==10000:
		print "Correct number of SNPs and participants"
	else:
		passTests=False;
		print "Error loading the genotype data!"

	if passTests:
		print "PASSES!!!"
	else:
		print "NOT PASS!!"

	print "\n\n\n"

##
##Some sanity checks for MUMAT
##
def testMUMAT():
	testfile1="TestCases/pop";	
	[y,sFil]=getData(testfile1);

	print "First we test the MU_Mat class!";
	try:
		mm=MU_Mat(sFil,[]);
	except:
		print "Trouble with MU_Mat construction!"
		return;
	print "A MU_Mat instance has been created!";
	try:
		mm.prod(y);
	except NotImplementedError:
		print "Next function passed!"
	else:
		print "A method in the MU_Mat class is not correct!";

	try:
		mm.score(y);
	except NotImplementedError:
		print "Next function passed!"
	else:
		print "A method (score) in the MU_Mat class is not correct!";
	


	try:
		mm.neighDist(y,1.0);
	except NotImplementedError:
		print "Next function passed!"
	else:
		print "A method in the MU_Mat class is not correct!";
	
	try:
		mm.sens(10);
	except NotImplementedError:
		print "Next function passed!"
	else:
		print "A method in the MU_Mat class is not correct!";


	try:
		mm.maxMU();
	except NotImplementedError:
		print "Next function passed!"
	else:
		print "A method in the MU_Mat class is not correct!";


	if mm.normY(y)!=[1.0,-1.0]:
		print "normY is broken in MU_Mat class";
	
	known=[1,2,3];

	try:
		snp=mm.snp_Names([1,2,3]);
	except:
		print "Error thrown by MU_Mat on snp_Names";
	snpReal=["null_1","null_2","null_3"];
	if snpReal!=snp:
		print "Wrond SNP names returned in MU_Mat";
	else:
		print "Correct SNP names!";

	try:
		ind=mm.snp_index(snp)
	except:
		print "Error thrown by MU_Mat in snp_index";
	if len([1 for i in range(0,len(ind)) if ind[i]==known[i]])!=len(ind):
		print "Incorrect indices returned!"
	else:
		print "snp_index returned the correct indices!";

	print "\n\n\n\n"
	print "Move on to test MU_Mem (assumes can do everythig in memory)!"

	try:
		mm=MU_Mem(sFil,[]);
	except:
		print "Error making MU_Mem object!!";
		return;
	else:
		print "Successfully created MU_Mem!!";

	if not isinstance(mm,MU_Mat):
		print "MU_Mem no longer extends MU_Mat!"
		print "Fail!"
		return;

	print "MU_Mem still extends MU_Mat!"

	if mm.BED!=sFil:
		print "Incorrect BED file!!";
		print "Fail!";
		return;
	if not isinstance(mm.X,np.ndarray):
		print "mm.X is not an array!"
	
	shp=np.shape(mm.X);
	print shp;
	if shp!=(10000,10000):
		print "X is wrong size!";
		return;
	print "X is correct!"

	print "Check maxMU"
	for l in range(0,20):
		MU=[[rand.uniform(-1,1) for j in range(0,20)] for i in range(0,100)]
		mm.MU=np.asarray(MU);
		mxMU=mm.maxMU();
		for i in range(0,100):
			if mxMU[i]!=max([abs(j) for j in MU[i]]):
				print "maxMU failed!"
				return;
	print "maxMU seems to pass!"

	print "Check Prod";

	for l in range(0,20):
		MU=[[rand.uniform(-1,1) for j in range(0,20)] for i in range(0,100)]
		y=[rand.uniform(0,1) for j in range(0,20)]
		mm.MU=np.asarray(MU);
		mxMU=mm.prod(y);
		for i in range(0,100):
			rl=sum([y[j]*MU[i][j] for j in range(0,20)]);
			if abs(mxMU[i]-rl)>.0001:
				print "prod failed!"
				return;
	print "prod seems to pass!"

	print "Check sens"


	for l in range(0,20):
		for mret in [i for i in range(0,10)]:
			mm.MU=np.asarray(MU);
			sen=mm.sens(mret);
			tru=0.0;
			for i in range(0,20):
				mu=[abs(MU[j][i]) for j in range(0,100)]
				mu=sorted(mu,reverse=True);
				tru=max(tru,sum(mu[:mret]))
			if sen!=tru:
				print "sens failed!"
				return;
	print "sens seems to pass!"


	print "Finally, check neighbor distance!";

	for l in range(0,20):
			MU=[[rand.uniform(-1,1) for j in range(0,100)]]
			mu=MU[0];
			mm.MU=MU;

			top=max(-sum([i for i in mu if i<0]),sum([i for i in mu if i>0]))
			y=[1.0 for i in range(0,100)]
			if not -100==mm.neighDist(y,top+1.0)[0]:
				print "Error in neigh dist!"
				return;	
	
			MU=[[rand.uniform(0,1) for j in range(0,100)]]
			mu=MU[0];
			mm.MU=MU;
			mu2=sorted(mu);
			val=min(mu)/2.0;
			for k in range(0,100):
				if (100-k)!=mm.neighDist(y,val)[0]:
					print "Error with neigh!";
					print k;
					return;
				val=val+mu2[k];	
	
			MU=[[rand.uniform(-1,1) for j in range(0,100)]]
			for k in range(0,20):
				y=[rand.randint(0,1) for i in range(0,100)];
				a=rand.randint(0,99);	
				x1=mm.neighDist(y,top/2.0)[0];
				y[a]=1-y[a];
				x2=mm.neighDist(y,top/2.0)[0];
				if abs(x1-x2)>1:
					print "Error in neigh dist!"
					return;
			

			
	print "neighdist appears ok!"

	print "MU_Mat seems ok!"
	print "\n\n\n"


##
##Tests MU_STRAT!
##
def TestMU_STRAT():
	
	print "Test MU_Strat!"
	testfile1="TestCases/pop";	
	[y,sFil]=getData(testfile1);
	sFil=sFil[:900,:1000]
	y=y[:900]

	mm=MU_STRAT(sFil,5);

	try:
		mm=MU_STRAT(sFil,5);
	except:
		print "Error creating MU_STRAT!"
		return;
	else:
		print "Created MU_STRAT!"

	if isinstance(mm,MU_Mem):
		print "Is a MU_Mem";
	else:
		print "Is not a MU_Mem!"
		return;


	print "Check Uk";
	try:
		Uk=mm.Uk
	except:
		print "Uk not generated!"
		return;
	else:
		print "Uk is generated!"

	a=np.shape(Uk)
	if a[0]!=900 and a[1]!=5:
		print "Error in dimensions of Uk!"
		return;
	print "Dimensions correct!"

	I=np.dot(Uk.T,Uk);
	for i in range(0,5):
		for j in range(0,5):
			if i!=j and abs(I[i][j])>.0001:
				print "Error in Uk!"
				print [i,j];
				return;
			elif i==j and abs(I[i][j]-1)>.0001:
				print "Error in Uk!"
				print i;
				return
	print "Uk seems to conists of k unit, orthogonal vectors!"

	mm.X=np.diag([10.0-i for i in range(0,10)]);

	for k in range(0,10):
		mm.calcMU(k);
		Uk=mm.Uk;
		for i in range(0,10):
			for j in range(0,k):
				if i!=j and abs(Uk[i][j])>.0001:
					print "Uk fails on diag!"
					return;
	print "Passes on diag!"
		

	print "Check normY";

	mm.X=np.diag([10.0-i for i in range(0,10)]);

	for k in range(1,10):
		mm.calcMU(k);
		for l in range(0,20):
			y=[rand.uniform(0,1) for i in range(0,10)];
			[bot,val]=mm.normY(y);
			mn=sum(y)/float(len(y));
			y=[i-mn for i in y];
			y=y[k:];
			val=sum([i**2 for i in y])
			if abs(bot**2-val)>.001:
				print "Error in normY!";
				print mn; 
				print bot**2;
				return;

	print "normY seems ok"

	print "Finally do some sanity checks on MU";

	for l in range(0,20):
		X=[[rand.uniform(-1,1) for i in range(0,50)] for j in range(0,100)]	
		
		sm=[sum(x)/50.0 for x in X]

		X=[[a-sm[i] for a in X[i]] for i in range(0,100)];
		X=np.asarray(X).T;
		mm.X=X;
		mm.calcMU(k);
		MU=mm.MU;
		Uk=mm.Uk;

		if max([abs(sum(mu)) for mu in MU])>.001:
			print "Error in calculating MU!"
			return;
	
		dt=np.dot(MU,Uk);
		if np.max(dt)>.001:
			print "Error in calculating MU!"
			return;
		if abs(max([sum([i**2 for i in mu]) for mu in MU])-1)>.001:
			print "Error in calculating MU!"
			return;	

		for i in range(0,100):
			x=[X[j][i] for j in range(0,50)];
			mu=MU[i]
			val=np.dot(x,mu);
			[bot,some,c]=mm.normY(x);	
			if abs(bot-val)>.001:
				print "Error in MU!"
				print bot;
				print val;
				return;	
	print "MU seems ok!"

		
	print "So MU_STRAT seems ok!!"



##
##Also want to test DP_UTIl!
##
def test_DP_UTIL():
	print "NEED TO IMPLEMENT";
	print "Compared results to old PrivSTRAT and seemed consistant";




##
##Test MU_LMM
##
def TestMU_LMM():	
	print "Test MU_LMM!"
	testfile1="../../GWAS/cleaned";	
	[y,sFil]=getData(testfile1);

	mm=MU_LMM(sFil,[1,-1.0]);

	try:
		mm=MU_LMM(sFil,[1,-1.0]);
	except:
		print "Error creating MU_LMM!"
		return;
	else:
		print "Created MU_LMM!"

	if isinstance(mm,MU_Mem):
		print "Is a MU_Mem";
	else:
		print "Is not a MU_Mem!"
		return;

	#mm=MU_LMM(sFil,[5,-1.0]);


	print "Do some sanity checks on MU";

	a=max(np.abs(np.sum(mm.MU,axis=1)))

	b=max(np.sum(np.abs(mm.MU),axis=1))

	if a/b>.001:
		print "Error with MU!"

	"""
	for l in range(0,20):
		X=[[rand.uniform(-1,1) for i in range(0,50)] for j in range(0,100)]	
		
		sm=[sum(x)/50.0 for x in X]

		X=[[a-sm[i] for a in X[i]] for i in range(0,100)];
		X=np.asarray(X).T;
		mm.X=X;
		mm.calcMU(k);
		MU=mm.MU;
		Uk=mm.Uk;

		if max([abs(sum(mu)) for mu in MU])>.001:
			print "Error in calculating MU!"
			return;
	
		dt=np.dot(MU,Uk);
		if np.max(dt)>.001:
			print "Error in calculating MU!"
			return;
		if abs(max([sum([i**2 for i in mu]) for mu in MU])-1)>.001:
			print "Error in calculating MU!"
			return;	

		for i in range(0,100):
			x=[X[j][i] for j in range(0,50)];
			mu=MU[i]
			val=np.dot(x,mu);
			[bot,some]=mm.normY(x);	
			if abs(bot-val)>.001:
				print "Error in MU!"
				print bot;
				print val;
				return;	
	print "MU seems ok!"
	"""
		
	print "So MU_LMM seems ok!!"




	

if __name__=="__main__":
	print "Let's get down to business!";
	#testloadfile();
	testMUMAT()
	TestMU_STRAT();
	#TestMU_LMM();
	test_DP_UTIL();

