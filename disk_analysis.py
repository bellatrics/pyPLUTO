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

class Force_Ana(objet):
    def Gravity(self,Data):
        [r2d, z2d] = np.meshgrid(Data.x1,Data.x2)
        r2d=r2d.T
        z2d=z2d.T
        Tool = pp.Tools()
        Grav_force_dict = {}
        Grav_force_dict['G_r'] = 1.0/(r2d*r2d)
        Grav_force_dict['G_z'] = np.zeros(r2d.shape)

        return Grav_force_dict

    def Pressure(self,Data,phi=10):
        Tool = pp.Tools()
        Prgrad = Tool.Grad(Data.pr[:,:,phi],Data.x1,Data.x2,Data.dx1,Data.dx2,polar=True)
        Press_force_dict ={}
        Press_force_dict['Fp_r'] = -1.0*(Prgrad[:,:,0]/Data.rho[:,:,phi])
        Press_force_dict['Fp_th'] = -1.0*(Prgrad[:,:,1]/Data.rho[:,:,phi])

        return Press_force_dict

    def Centrifugal(self,Data,phi=10):
        [r2d,z2d] = np.meshgrid(Data.x1,Data.x2)
        r2d=r2d.T
        z2d=z2d.T
        Centri_force_dict={}
        Centri_force_dict['Fcf_r'] = (Data.v3[:,:,10]*Data.v3[:,:,10])/r2d
        Centri_force_dict['Fcf_z'] = np.zeros(r2d.shape)
	
        return Centri_force_dict


        
    
    
    
