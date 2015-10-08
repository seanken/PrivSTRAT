###############################
##Some useful methods
###############################
from numpy.linalg import svd
import math;
import numpy as np;
#import matplotlib.pyplot as plt;
from numpy import random as rand;

##
##Gen Laplacian noise
##
def Lap(mu,lamb):
	return rand.laplace(mu,lamb);	

##
##normalize X
##
def norm(X):
	return X;
	n=len(X);
	m=len(X[0]);
	NX=np.asmatrix(X).T;
	maf=np.mean(X,axis=0);
	bot=[math.sqrt(mi*(1-mi/2.0)) for mi in maf];
	NX=X-maf;

	NX=NX/np.asarray(bot);

	return NX;
##
##Samples from a normal with mean mu, variance var
##Returns n samples
##
def N(mu,var,n=1):
	nrm=rand.normal(mu,math.sqrt(var),n)[:];
	return [nr for nr in nrm]

##
##Finds the number of elements in both a and b
##
def inter(a,b):
	ln=float(len(a));
	return len([i for i in a if i in b])/ln;

"""
##
##Plots stuff
##
def PlotIt(x,ys,ytitle="",xtitle="",txtSize=16,cols=[],saveName="",bnds=False):
	fig=plt.figure();
	ax=fig.add_subplot(1,1,1)
	i=0;
	for y in ys:
		if len(cols)==0:	
			ax.plot(x,y);
		else:
			ax.plot(x,y,cols[i]);
			print cols[i];
			i=i+1;
	ax.spines['top'].set_visible(False)

	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('none');

	ax.yaxis.set_ticks_position('none');

	plt.xlabel(xtitle,fontsize=20);
	plt.ylabel(ytitle,fontsize=20);	
	if bnds:
		plt.ylim(0,1.0);
	if not bnds:
		plt.ylim(0,int(max([max(y) for y in ys])));
	if len(saveName)==0:
		plt.show();
	else:
		fig.savefig(saveName);

"""	
