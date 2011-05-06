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


def disk_reso(Data,ul=None):
    print " This pro requires to have unit length of simumation."
    print "By Default its taken to be 1 AU"
    if ul==None : ul=1.0
    
    Delr = Data.dx1.max()*ul
    Deltheta = Data.x1[Data.n1-1]*Data.dx2.min()*ul
    DelPhi = Data.x1[Data.n1-1]*np.sin(Data.x2[0.5*Data.n2])*Data.dx3.min()*ul
    Rat1 = Delr/Delr
    Rat2 = Deltheta/Delr
    Rat3 = DelPhi/Delr
    print '-----------ANALYSIS RESULTS----------------'
    print 'Delta R = %f in the units of scale length ul'%(Delr)
    print 'Delta theta = %f in the units of scale length ul'%(Deltheta)
    print 'Delta phi = %f in the units of scale length ul'%(DelPhi)
    print 'THUS the ratio is :%f:%f:%f'%(Rat1,Rat2,Rat3)
    print '--------------------------------------------'



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

        
    
    
    
