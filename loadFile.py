######################
##Loads Plink file
######################
import numpy as np;
from pysnptools.snpreader import Bed;
from pysnptools.snpreader import Pheno;
from os.path import isfile;

##
##gets y, the phenotype, X, the normalized matrix, sFil,
##the snpreader, and Q the covariate matrix
##
def getData(filename):
	mph=3;
	sFil=Bed(filename);
	yFil=Pheno(filename+".fam");
	
	y=yFil.read().val[:,mph];
	y=[i-1 for i in y]
	return [y,sFil];
	


