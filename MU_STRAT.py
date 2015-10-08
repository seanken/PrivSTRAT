###########
##Implements the MU matrix class for
##PrivSTRAT (see manuscript)
###################
from MU_Mat import MU_Mat;
from MU_Mat import MU_Mem;
from sklearn.decomposition import TruncatedSVD as svdt;
import math;
import numpy as np;

class MU_STRAT(MU_Mem):
    ##
    ##Implementation of calc MU matrix for PrivSTRAT
    ##
    def calcMU(self,k):
        self.k=k
	uk_temp=svdt(n_components=k);
        u=np.asarray(uk_temp.fit_transform(self.X));
        bot=np.sum(u**2,axis=0);
        bot=[math.sqrt(i) for i in bot]
        self.Uk=u/bot;
    
        MU_temp=self.X-np.dot(self.Uk,np.dot((self.Uk).T,self.X));
    
        self.MU=MU_temp.T;
    
        sd=np.sum(self.MU**2,axis=1);
        sd=[math.sqrt(s) for s in sd];
        self.MU=self.MU/np.asarray(sd)[:,np.newaxis];
    

    ##
    ##Returns list of scores
    ##
    def score(self,y):
        MUy=prod(y);
        bot=1.0;
        return [my**2/bot**2 for my in MUy];

    ##
    ##normalizes y, by mean centering and projecting out U_k,
    ##then returns the resulting vectors length
    ##
    def normY(self,y):
        n=len(y);
        mn=sum(y)/float(len(y));
        y_st=[i-mn for i in y];
        y_st=np.asarray(y_st)-np.dot(self.Uk,np.dot(self.Uk.T,y_st));
        return [math.sqrt(sum([i**2 for i in y_st])),math.sqrt(float(n-1)/float(n)),math.sqrt(n-self.k-1)]
