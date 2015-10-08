###########
##Implements the MU matrix class for
##PrivLMM (see manuscript)
###################
import numpy as np;
from MU_Mat import MU_Mat;
from MU_Mat import MU_Mem;
import math;
from numpy.linalg import inv;

class MU_LMM(MU_Mem):
	##
	##Implementation of calc MU matrix for PrivLMM
	##
	def calcMU(self,[se2,sg2]):
		n=len(X);
		m=len(X[0])
		print "Taking inverse!"
		Kinv=inv(np.eye(n)*se2+sg2/float(m)*np.dot(self.X,(self.X).T))
   		print "Calc Top"
		self.MU=np.dot(X.T,Kinv);
		print "Calc Bottom"
		bot=np.asarray([np.dot(self.X[:,i],self.MU[i]) for i in range(0,m)])
		print "Divde through!!"
		self.MU=self.MU/bot[:,np.newaxis]	 
		print "Done!";

	##
	##Calculates variance using ML estimate, 
	##If epsilon>0 does so in epsilon-DP way (see Supplement to manuscript).
	##
	def calcVar(BED,num=1.0,epsilon=-1.0):
		


