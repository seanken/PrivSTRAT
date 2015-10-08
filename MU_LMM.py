###########
##Implements the MU matrix class for
##PrivLMM (see manuscript)
###################
import numpy as np;
from MU_Mat import MU_Mat;
from MU_Mat import MU_Mem;
import math;
from numpy.linalg import inv
import math;
from numpy.random import permutation as perm
import os;
from pysnptools.util import intersect_apply;
from pysnptools.snpreader import Bed;
from pysnptools.snpreader import Pheno;
from os.path import isfile;
from numpy.random import laplace as Lap;
from fastlmm.inference import LMM;


class MU_LMM(MU_Mem):
	##
	##Implementation of calc MU matrix for PrivLMM
	##
	def calcMU(self,par):
		num=par[0];
		epsilon=par[1];
		if len(par)>2:
			self.VarCalc=par[3];
		else:
			self.VarCalc=FastLMM();
		n=len(self.X);
		m=len(self.X[0])
		print "Calculate variance"
		[se2,sg2]=self.estVar(num,epsilon)
		self.se2=se2;
		self.sg2=sg2;
		print "Taking inverse!"
		Kinv=inv(np.eye(n)*se2+sg2/float(m)*np.dot(self.X,(self.X).T))
   		print "Calc Top"
		self.MU=np.dot(self.X.T,Kinv);
		print "Calc Bottom"
		bot=np.asarray([np.dot(self.X[:,i],self.MU[i]) for i in range(0,m)])
		print "Divide through!!"
		self.MU=self.MU/bot[:,np.newaxis]
		print "Center!"
		mn=np.mean(self.MU,axis=1);
		print len(mn); 
		self.MU=self.MU-mn[:,np.newaxis];
		print "Done!";


		




	##
	##estimates the variance of y
	##
	def estVarY(self,y,epsilon):
		vr=np.var(y);
		n=len(y);
		return vr+Lap(0,3/float(epsilon*n))

	##
	##Divide data!
	##
	def divideData(self,filename,num=5,mph=3,delet=True):
		print "Estimating heritability using "+str(num)+" components"
		direct="TEMP"
		sFil=Bed(filename);
		yFil=Pheno(filename+".fam");
		n=sFil.iid_count	
		reOrd=perm(n);
		yFil=yFil[reOrd,:];
		sFil=sFil[reOrd,:];

		y=yFil.read().val[:,3];

		div=[int(math.ceil( i*n/float(num) )) for i in range(0,num+1)];
		
		varEsts=[];

		for i in range(0,num):
			print "For component "+str(i);
			sFilTemp=self.BED[div[i]:div[i+1],:];
			Xtemp=sFilTemp.read().standardize().val;
			ytemp=y[div[i]:div[i+1]];

			varEsts.append(self.VarCalc.RealVar(ytemp,Xtemp));
		
		return varEsts;

	##
	##Estimate Variance in \epsilon-DP way!
	##
	def estVar(self,num,epsilon):
		filename=self.BED.filename;
		y=Pheno(filename+".fam").read().val[:,3];
		varEsts=self.divideData(filename,num=num);
		if epsilon<0:
			return varEsts[0];
		e1=.1*epsilon;
		e2=.45*epsilon;
		e3=.45*epsilon;
		vary=self.estVarY(y,e1);
		se2=sum([v[1] for v in varEsts])/float(num)+Lap(0.0,vary/(e2*float(num)));
		if se2<0:
			se2=0;
		if se2>vary:
			se2=vary;
		sg2=sum([v[0] for v in varEsts])/float(num)+Lap(0.0,vary/(e3*float(num)));

		if sg2<0:
			sg2=0;
		if sg2>vary:
			sg2=vary;

		return [sg2,se2];


##
##Class used to estimate variance
##
class VarEstimator:
	def realVar(self,y,X):
		raise NotImplementedError("Not implemented!");

from fastlmm.inference import LMM;

##
##Uases FASTLMM
##
class FastLMM():
	##
	##uses REML to estimate variance components
	##
	def RealVar(self,y,X):
		lmmg=LMM()
		m=np.shape(X)[1];
		n=len(y);
		lmmg.setG(X/math.sqrt(m))
		lmmg.sety(y);
		lmmg.setX(np.ones([n,1]))
		dct=lmmg.findH2();
		h2=dct['h2'];
		s2=dct['sigma2'];
		sg2=h2*s2;
		se2=s2-sg2;
		return [se2,sg2];


