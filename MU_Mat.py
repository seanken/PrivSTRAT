###################
##Implements the methods related to the MU matrix
#####################
import pysnptools;
import numpy as np;
import scipy as sp;



##########################
##This class represents
##the MU matrix required for our method
##note that this version is abstract, most
##of the methods are not implemented
##This is so that, later on, we can improve our code so it
##can work on genotype data that can not be read into
##memory due to size constraints
#########################

class MU_Mat:
	##
	##Initialize
	##
	def __init__(self,BED,par,binSize=1):
		self.BED=BED;
		self.makeBins();
	
	def makeBins(self,binSize=1):
		self.binSize=binSize;
		self.Pos=self.BED.pos[:,0:3:2]
		self.Pos=[[p[0],int(p[1]/binSize)] for p in self.Pos];
		if binSize>1:
			I=sorted([i for i in range(0,self.BED.sid_count)],key=lambda i:self.Pos[i])
			self.BED=self.BED[:,I];		
			self.Pos=sorted(self.Pos);
	

	
	##
	##Returns the product of MU with y
	##where y is a 
	##
	def prod(self,y):
		raise NotImplementedError("Need to specify how to calculate product!");
	
	##
	##Returns a list the \chi^2 stat for all SNPs
	##Differs between methods so not implemented here
	##
	def score(self,y):
		raise NotImplementedError("Need to specify how to calculate Wald Score!");
	##
	##Returns [b_i(c) i a SNP]
	##See manuscript for details
	##
	def neighDist(self,y,c):
		raise NotImplementedError("Need to specify how to calculate the neighbor distance!");

	##
	##Given mret, returns the associated sensitivity value, \Delta
	##(See manuscript for details)
	##
	def sens(self,mret):
		raise NotImplementedError("Need to specify how to calculate sensitivity!");

    	##
    	##Returns a list of max_j|mu_{ij}| for all i
    	##
   	def maxMU(self):
        	raise NotImplementedError("Need to specify how to calculate product!");


    	##
   	##returns a normalization constant, based on y
    	##As well as sensitivity
    	##
    	def normY(self,y):
        	return [1.0,-1.0]


    	##
    	##Get the names of SNPs at specified indices
    	##
    	def snp_Names(self,ind,binned=False):
		if binned:
			return [self.Pos[i] for i in ind]
	     	return [self.BED.sid[i] for i in ind];

    	##
    	##Gets indices of SNPs with specified names
    	##
    	def snp_index(self,snps):
        	return self.BED.sid_to_index(snps);


##
##An extension of MU_Mat that assumes that MU and X can fit in memory
##
class MU_Mem(MU_Mat):	
	##
	##Initialize, where BED is a 
	##MU used if already calculated MU
	##SNPs specified if only looking at certain SNPs
	##
	def __init__(self,BED,par,binSize=1):
		self.BED=BED
		self.makeBins(binSize);	

		self.X=BED.read().standardize().val##normalized genotype data
        	self.sensit=-1.0;	
		self.MUn=[];
		self.MUp=[];
		self.y=[];
		self.calcMU(par);


	##
	##Used to calculate the actual array MU
	##
	def calcMU(self,par):
		self.MU=[];

	##
	##Returns the product of MU with y
	##where y is a 
	##
	def prod(self,y):
		return np.dot(self.MU,y);		
	

	##
	##Returns [b_i(c) i a SNP]
	##See manuscript for details
	##If reuse=True means that neighbor distance
	##will be calculated numerous times for same y and MU but different c
	##so reuse some calculations (though comes at the cost of more memory);
	##
	def neighDist(self,y,c,reuse=False):
		bnd=c;
		m=len(self.MU);
		n=len(self.MU[0]);
		sc2=np.dot(self.MU,y);#.flatten();
		sc=np.abs(sc2);
		if reuse and len(self.MUn)>0:
			MUn=self.MUn;
			MUp=self.MUp;
		else:
			MUz=self.MU*(1-2*np.asarray(y));
			MUp=np.maximum(MUz,0);
			MUn=np.minimum(MUz,0);
			MUp=np.sort(-MUp,axis=1);
			MUn=np.sort(MUn,axis=1);
			MUp=-np.cumsum(MUp,axis=1);
			MUn=np.cumsum(MUn,axis=1);
			if reuse:
				self.MUn=MUn;
				self.MUp=MUp;
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
	##Given mret, returns the associated sensitivity value, \Delta
	##(See manuscript for details)
	##
	def sens(self,mret):
		#if self.sensit>0:
		#	return self.sensit;
		MU2=np.absolute((self.MU).T)				
		MU2=np.sort(-MU2)[:,:mret];
		sens=min(np.sum(MU2,axis=1));
		sens=-sens;
		self.sensit=sens;
		return sens;
	

	##
    	##Returns a list of max_j|mu_{ij}| for all i
    	##
    	def maxMU(self):
        	return np.max(np.absolute(self.MU),axis=1)


