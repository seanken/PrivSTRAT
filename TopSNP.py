###################
##Implements methods to pick top SNPs
#####################
#import matplotlib.pyplot as plt;
import numpy as np;
import math;
import random as rand;
from numpy.random import multinomial as mult;
import Useful as us;

##
##Returns mret distinct elements in the range
##0 to len(sc)-1, inclusive
##Picks i with prob prop to exp(-sc[i])
##
def expPick(sc,mret):
	m=len(sc);
	fnd=[0 for i in range(0,mret)]
	scores=[s for s in sc]
	for i in range(0,mret):
		ms=max(scores);
		scores=[s-ms for s in scores];
		prb=[math.exp(s) for s in scores];
		sm=sum(prb);
		prb=[s/sm for s in prb];
		while sum(prb)>1:
			sm=sum(prb);
			prb=[s/sm for s in prb];

		mlt=mult(1,prb);

		k=min([k for k in range(0,m) if mlt[k]>0]);
		
		fnd[i]=k;
		scores[k]=-float("inf");
	return fnd;



##
##An optimized version of the method to find the neighbor distance
##
def fastNeigh(bnd,MU,y):
	m=len(MU);
	n=len(MU[0]);
	sc2=np.dot(MU,y);#.flatten();
	sc=np.abs(sc2);
	MUz=MU*(1-2*np.asarray(y));
	MUp=np.maximum(MUz,0);
	MUn=np.minimum(MUz,0);
	MUp=np.sort(-MUp,axis=1);
	MUn=np.sort(MUn,axis=1);
	MUp=-np.cumsum(MUp,axis=1);
	MUn=np.cumsum(MUn,axis=1);
	cntp=[np.searchsorted(MUp[i],bnd-sc2[i]) for i in range(0,m)];
	cntn=[np.searchsorted(-MUn[i],bnd+sc2[i]) for i in range(0,m)];
	neigh=-np.minimum(cntp,cntn);
	neigh=[i for i in neigh];
	I=[i for i in range(0,m) if sc[i]>bnd];
	for i in I:
		mu=[]
		if sc2[i]>bnd:
			mu=MUn[i];
			j=np.searchsorted(-mu,-bnd+sc2[i]);
			neigh[i]=j+1;
		if sc2[i]<-bnd:
			mu=MUp[i];
			j=np.searchsorted(mu,-bnd-sc2[i]);
			neigh[i]=j+1
	return neigh;
	
	

	




##
##Picks top scoring SNPs using neighbor method
##
def PickTopNeigh(y,MU,mret,epsilon,neighDist=[],randBND=True,bnd=-1):
	n=len(y);
	m=len(neighDist);
	ep1=.1*epsilon;
	ep2=.9*epsilon;
	if len(neighDist)==0:
		sc=np.dot(MU,y);#[abs(sum([y[i]*mu[i] for i in range(0,n)])) for mu in MU];
		sc=[abs(s) for s in sc]

		if bnd<0:
			bnd=sum(sorted(sc,reverse=True)[mret-1:mret+1])/2.0;
			if randBND:
				bnd=bnd+us.Lap(0,np.max(np.abs(MU))/ep1);
				bnd=abs(bnd);
		print "Calculating Distance"
		neighDist=fastNeigh(bnd,MU,y)
		m=len(MU);

	sc=[nei*ep2/(2.0*mret) for nei in neighDist];

	return expPick(sc,mret);


##
##Picks top scoring SNPs using Score method
##
def PickTopScore(y,MU,mret,epsilon,sc,sens):
	n=len(y);
	m=len(sc);
	ms=max(sc);
	sc=[s-ms for s in sc];

	sc2=[s*epsilon/(2.0*sens) for s in sc];
	
	return expPick(sc2,mret);




##
##Picks top scoring SNPs using Laplacian method
##
def PickTopLap(y,MU,mret,epsilon=1.0,sc=[],sens=-1):
	n=len(y);
	m=len(sc);
	deltQ=sens;

	scDP=[s+us.Lap(0,2*deltQ/epsilon) for s in sc];

	return sorted([i for i in range(0,m)],key=lambda i:-scDP[i])[:mret];


##
##epsilon-DP estimate of number
##Uses neigh dist
##
def estNum(neigh,epsilon):
	m=len(neigh);
	if epsilon<0:
		return len([i for i in neigh if i>0]);
	sc=[0.0 for i in range(0,m+1)];
	for i in range(0,m):
		if neigh[i]>0:
			sc[i]=-epsilon*neigh[i]/2.0;
		else:
			sc[i]=epsilon*(neigh[i]-1)/2.0;
	for i in range(0,m+1):
		sc[i]=math.exp(sc[i]);
	sm=sum(sc);
	sc=[i/sm for i in sc]
	multi=mult(1,sc);
	ret=min([i for i in range(0,m+1) if multi[i]>0]);
	if ret==m:
		return len([i for i in neigh if i>0]);
	val=neigh[ret];
	
	v1=len([i for i in neigh if i>val]);
		
	v2=len([i for i in neigh if i==val]);


	return rand.randint(v1,v1+v2-1);


##
##Calcs the sensitivity
##
def getSens(MU,mret):
	MU2=np.absolute(np.asarray(MU).T);
	MU2=np.sort(-MU2)[:,:mret];
	sens=min(np.sum(MU2,axis=1));
	sens=-sens;
	print "Calculated Sensitivity!"
	return sens;
