import pyPLUTO as pp
import PhyConst as phc
import os
import numpy as np
import asciidata as asc
from scipy import integrate


##############--DATA ANALYSIS CLASS--##########################
#
# This Class has all the functions for the analysing the Data.
# and getting fancy quantities from data.
# CALLED AFTER pyPLUTO.pload object is defined.
#
################################################################

def Analysis(filepath=None,info=None):
    Values = np.asarray(asc.open(filepath+'analysis.out'))
    [a,b]=Values.shape
    
    MyArr = []
    finfo = open(filepath+'analysis.info','r')
    for line in finfo.readlines():
        if line.find('column')>=0:
            MyArr.append(line.split())

    print len(MyArr)
    ana_dict= dict([('col'+ str(i+1),[Values[i],' '.join(MyArr[i][3:])]) for i in range(a-1)])
    
    
    
    return ana_dict

        
    
    
    
