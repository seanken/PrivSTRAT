###################
##method to return top SNPs
##in privacy preserving manner
#########################
from MU_Mat import MU_Mat;
from scipy.stats import chi2;
import numpy as np;
import scipy as sp;
from numpy.random import multinomial as mult;
import random as rand;
import math;
from numpy.random import laplace as Lap;




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
##Picks top scoring SNPs using neighbor method
##reuse indicates will calc neigh dist lost of times for same y, 
##different bnds, so saves some calculations to speed it up
##(mostly for making figures faster)
##
def PickTopNeigh(y,MU,mret,epsilon,reuse=False):
    n=len(y);
 
    ep1=.1*epsilon;
    ep2=.9*epsilon;

    sc=MU.prod(y);
    sc=[abs(s) for s in sc]
                
    bnd=sum(sorted(sc,reverse=True)[mret-1:mret+1])/2.0;
    
    bnd=bnd+Lap(0,np.max(np.abs(MU.MU))/ep1);
    bnd=abs(bnd);
    print "Calculating Distance"
    neighDist=MU.neighDist(y,bnd,reuse=reuse);
    
    
    sc=[nei*ep2/(2.0*mret) for nei in neighDist];
        
    index_Ret=expPick(sc,mret);

    return MU.snp_Names(index_Ret);

##
##Picks top scoring SNPs using Score method
##
def PickTopScore(y,MU,mret,epsilon):
    n=len(y);
    sc=MU.prod(y);
    sc=[abs(s) for s in sc]
    m=len(sc);
    sens=MU.sens(mret);
    ms=max(sc);
    sc=[s-ms for s in sc];
        
    sc2=[s*epsilon/(2.0*sens) for s in sc];
        
    index_Ret=expPick(sc2,mret);

    return MU.snp_Names(index_Ret);

##
##Picks top scoring SNPs using Laplacian method
##
def PickTopNoise(y,MU,mret,epsilon):
    n=len(y);
    sc=MU.prod(y);
    sc=[abs(s) for s in sc]
    m=len(sc);
    sens=MU.sens(mret);

    if epsilon<0:
	scDP=sc;
    else:       
    	scDP=[s+Lap(0,2*sens/epsilon) for s in sc];
        
    index_Ret=sorted([i for i in range(0,m)],key=lambda i:-scDP[i])[:mret];

    return MU.snp_Names(index_Ret);




##
##The interface for the outside world
##algor can be either noise, score or neighbor
##
def PickTop(y,MU,mret,epsilon,algor="noise",reuse=False):
    if algor=="noise":
        return PickTopNoise(y,MU,mret,epsilon);
    elif algor=="neighbor":
        return PickTopNeigh(y,MU,mret,epsilon,reuse=reuse);
    elif algor=="score":
        return PickTopScore(y,MU,mret,epsilon);
    else:
        raise ValueError("The algorithm "+algor+" is not recognized as a valid choice.");






##################
#Now deals with Wald Test
#########################
##
##Returns the estimate of the wald test for chosen SNPs using Laplacian mechanism
##
def WaldTest(y,MU,epsilon,snps,forFigs=False):
    if len(snps)==0:
            I=[i for i in range(0,np.shape(MU.MU)[0])]
	    eps=epsilon
    else:
	    I=MU.snp_index(snps);
            eps=epsilon/float(len(snps));

    if forFigs:
	eps=epsilon;

    sc=MU.prod(y);
    [nm,sen]=MU.normY(y);
    mxMU=MU.maxMU();
    
    bot=nm**2;

    if sen>0 and epsilon>0:
        eps=.5*eps;
        bot=(nm+Lap(0.0,2.0*sen/epsilon))**2;
    if forFigs:
	eps=epsilon/2.0;
	bot=[(nm+Lap(0,sen/eps))**2 for i in I]

    if eps>0:
    	sc=[sc[i]+Lap(0.0,mxMU[i]/eps) for i in I];


    if forFigs:
	return [sc[i]**2/bot[i] for i in range(0,len(sc))]
    return [s**2/bot for s in sc];








############################
##Finally deals with estimating number of significant SNPs
#########################
##
##epsilon-DP estimate of number
##Uses neigh dist
##
def estNum(MU,y,pval,epsilon):
    bnd_sc=math.sqrt(chi2.ppf((1.0-pval),df=1));

    ret=MU.normY(y);
    nm=ret[0];
    sen=ret[1];


    n=len(y);
    if sen>0:
        bnd_est=bnd_sc*abs(nm+Lap(0.0,2.0*sen/epsilon));
	if MU.k>0:
		bnd_est=bnd_est/math.sqrt(n-MU.k-1);
    else:
        bnd_est=bnd_sc*nm;

    
    neigh=MU.neighDist(y,bnd_est);
    m=len(neigh);
    if epsilon<0:
        return len([i for i in neigh if i>0]);
    sc=[0.0 for i in range(0,m+1)];
    for i in range(0,m):
        if neigh[i]>0:
            sc[i]=-epsilon*neigh[i]/4.0;
        else:
            sc[i]=epsilon*(neigh[i]-1)/4.0;
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



#####
##CI for wald
#########
##
##
##
def interLap(v,sens,p):
	r=-sens*math.log(1-p);
	a=v-r;
	b=v+r;
	if a<0 and b<=0:
		return [b**2,a**2];
	if a>=0 and b>0:
		return [a**2,b**2];
	return [0,max(a**2,b**2)];

##
##Differentially private interval estimate, p% CI
##
def CI(y,MU,p,epsilon,snps):
    I=MU.snp_index(snps);
    eps=epsilon/float(len(snps));
    q=math.sqrt(p);
    sc=MU.prod(y);
    [nm,sen]=MU.normY(y);
    mxMU=MU.maxMU();
    

    eps=.5*eps;
    botRoot=(nm+Lap(0.0,2.0*sen/epsilon));
    bot=botRoot**2
    sc=[sc[i]+Lap(0.0,mxMU[i]/eps) for i in I];

    yInter=interLap(botRoot,2.0*sen/epsilon,q);
    scInter=[interLap(sc[i],mxMU[I[i]]/eps,q) for i in range(0,len(I))];
    CIup=[scInter[i][1]/yInter[0] for i in range(0,len(I))]
    CIdown=[scInter[i][0]/yInter[1] for i in range(0,len(I))]
    waldEst=[s**2/bot for s in sc];
    return [CIdown,waldEst,CIup];


