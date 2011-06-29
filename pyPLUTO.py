import os
import sys
import struct
import numpy as np
import linecache
import scipy as S
import scipy.ndimage
from scipy import integrate
import scipy.interpolate
from matplotlib.pyplot import *

def curdir():
        curdir = os.getcwd()+'/'
        return curdir

def get_nstepstr(ns):
    nstepstr = str(ns)
    while len(nstepstr) < 4:
        nstepstr= '0'+nstepstr
    return nstepstr




##############--PLOAD CLASS--##########################
#
# This Class has all the routines loading the data from the
# binary files output from PLUTO Simulations. Assign an object
# when the data is loaded for some NSTEP
# Eg --
# import pyPLUTO as plp
# x = plp.pload(200) : Where 200 is data.200.dbl or NSTEP=200.
# Thus x is the pyPLUTO.pload object to play around.
#
################################################################


class pload(object):
    def get_varinfo(self,w_dir=None):
        if w_dir is None: w_dir =curdir()
        fname_v = w_dir+"dbl.out"
        f_var = open(fname_v)
        lnum_v = len(f_var.readlines())
        Var_info = linecache.getline(fname_v,1).split()
        fltype = Var_info[4]
        nvar = len(Var_info[6:])
        allvars=[]
        for i in Var_info[6:]:
            allvars.append(i)
        
            if i == 'b1s':
                nvar = nvar - 1
                del allvars[-1]
            
            if i == 'b2s':
                nvar = nvar - 1
                del allvars[-1]
            
                
        f_var.close()
        return {'fltype':fltype, 'nvar':nvar, 'allvars':allvars}

    def time_info(self,w_dir=None):
	    if w_dir is None: w_dir =curdir()
	    fname_v = w_dir+"dbl.out"
	    last_line = file(fname_v,"r").readlines()[-1].split()
	    nlast = int(last_line[0])
	    SimTime =  float(last_line[1])
	    Dt = float(last_line[2])
	    Nstep = int(last_line[3])
	    
	    print "------------TIME INFORMATION--------------"
	    print 'nlast =',nlast
	    print 'time  =',SimTime
	    print 'dt    =', Dt
	    print 'Nstep =',Nstep
	    print "-------------------------------------------"
	    
	    return {'nlast':nlast,'time':SimTime,'dt':Dt,'Nstep':Nstep}
    
	    
	    
    def geometry(self,w_dir=None):
        if w_dir is None : w_dir = curdir()
	fname_d = w_dir+"definitions.h"
        f_def = open(fname_d)
        lnum_d = len(f_def.readlines())
        Geo_info = linecache.getline(fname_d,4).split()
        f_def.close()
        print "GEOMETRY >> " + Geo_info[2]
        

    def grid(self,w_dir=None):
        if w_dir is None: w_dir=curdir()
        fname_g = w_dir+"grid.out"
        f_grid=open(fname_g)
        lnum_g = len(f_grid.readlines())
        n1 = linecache.getline(fname_g,1)
        n2 = linecache.getline(fname_g,int(n1)+2)
        n3 = linecache.getline(fname_g,int(n1)+int(n2)+3)
        x1=[]
        x2=[]
        x3=[]
        dx1=[]
        dx2=[]
        dx3=[]
        for i in range(2,int(n1)+2):
            A = linecache.getline(fname_g,i).split()
            x1.append(float(A[2]))
            dx1.append(float(A[4]))
        
        x1 = np.asarray(x1)
        dx1 = np.asarray(dx1)

        for j in range(3+int(n1),int(n1)+int(n2)+3):
            B = linecache.getline(fname_g,j).split()
            x2.append(float(B[2]))
            dx2.append(float(B[4]))
        
        x2 = np.asarray(x2)
        dx2 = np.asarray(dx2)

        for k in range(4+int(n1)+int(n2),lnum_g+1):
            C = linecache.getline(fname_g,k).split()
            x3.append(float(C[2]))
            dx3.append(float(C[4]))

        x3 = np.asarray(x3)
        dx3 = np.asarray(dx3)

        f_grid.close()
	
        grid_dict={'n1':int(n1),'n2':int(n2),'n3':int(n3),'x1':x1,'x2':x2,'x3':x3,'dx1':dx1,'dx2':dx2,'dx3':dx3}

	return grid_dict


    def data(self,ns,w_dir=None):
        if w_dir is None : w_dir=curdir()
	print "Working Dir : %s" % (w_dir)
        grid_dict = self.grid(w_dir)
        nstep = get_nstepstr(ns)
        varinf= self.get_varinfo(w_dir)
        data_dict={}
        n1 = grid_dict.get('n1')
        n2 = grid_dict.get('n2')
        n3 = grid_dict.get('n3')
        print "<DOMAIN> %d x %d x %d " % (n1,n2,n3)
        
        if varinf.get('fltype') == 'single_file':
            fname_data = w_dir+"data."+nstep+".dbl"
            f_data = open(fname_data,'rb')
            datout = f_data.read()
            D=struct.unpack("<"+str(len(datout)/8)+"d",datout)
            A = np.asarray(D)
            for i in range(varinf.get('nvar')):
                print "> Reading %s" % (varinf.get('allvars')[i])
                if varinf.get('allvars')[i] == varinf.get('allvars')[-1] :
			if n3 == 1:
				data_dict[(varinf.get('allvars')[i])]=A[-n2*n1:].reshape(n2,n1).transpose()
			else:
				data_dict[(varinf.get('allvars')[i])]=A[-n3*n2*n1:].reshape(n3,n2,n1).transpose()
                else :
			if n3 == 1:
				data_dict[(varinf.get('allvars')[i])]=A[i*n2*n1:(i+1)*n2*n1].reshape(n2,n1).transpose()
			else:
				data_dict[(varinf.get('allvars')[i])]=A[i*n3*n2*n1:(i+1)*n3*n2*n1].reshape(n3,n2,n1).transpose()
            
        else:
            fname_list = []
            f_list = []
            datout =[]
            Dind=[]


            for j in range(varinf.get('nvar')):
                fname_list.append(w_dir+ varinf.get('allvars')[j]+"."+nstep+".dbl")
                f_list.append(open(fname_list[j],'rb'))
                datout.append(f_list[j].read())
                Dind.append(struct.unpack("<"+str(len(datout[j])/8)+"d",datout[j]))

	    
	    A = np.asarray(Dind)
	    
	    for j in range(varinf.get('nvar')):
		    if n3 == 1:
			    data_dict[(varinf.get('allvars')[j])]=A[j].reshape(n2,n1).transpose()
		    else:
			    data_dict[(varinf.get('allvars')[j])]=A[j].reshape(n3,n2,n1).transpose()
		    

	
	return data_dict
                
    def __init__(self,ns,w_dir=None):
	    if w_dir is None : w_dir=curdir()
	    Grid_dictionary=self.grid(w_dir)
	    Data_dictionary=self.data(ns,w_dir)
	    for keys in Grid_dictionary:
		    object.__setattr__(self,keys,Grid_dictionary.get(keys))
	    for keys in Data_dictionary:
		    object.__setattr__(self,keys,Data_dictionary.get(keys))



##############--TOOLS CLASS--##########################
#
# This Class has all the functions doing basic mathematical
# operations to the vector or scalar fields.
# CALLED AFTER pyPLUTO.pload object is defined
#
################################################################


class Tools(object):
	def deriv(self,Y,X=None):
		n = len(Y)
		n2 = n-2
		if X==None : X = np.arange(n)
		Xarr = np.asarray(X,dtype='float')
		Yarr = np.asarray(Y,dtype='float')
		x12 = Xarr - np.roll(Xarr,-1)   #x1 - x2
		x01 = np.roll(Xarr,1) - Xarr    #x0 - x1
		x02 = np.roll(Xarr,1) - np.roll(Xarr,-1) #x0 - x2
		DfDx = np.roll(Yarr,1) * (x12 / (x01*x02)) + Yarr * (1./x12 - 1./x01) - np.roll(Yarr,-1) * (x01 / (x02 * x12))
		# Formulae for the first and last points:

		DfDx[0] = Yarr[0] * (x01[1]+x02[1])/(x01[1]*x02[1]) - Yarr[1] * x02[1]/(x01[1]*x12[1]) + Yarr[2] * x01[1]/(x02[1]*x12[1])
		DfDx[n-1] = -Yarr[n-3] * x12[n2]/(x01[n2]*x02[n2]) + Yarr[n-2]*x02[n2]/(x01[n2]*x12[n2]) - Yarr[n-1]*(x02[n2]+x12[n2])/(x02[n2]*x12[n2])

		return DfDx
	
	def Grad(self,phi,x1,x2,dx1,dx2,polar=False):
		(n1, n2) = phi.shape 
		grad_phi = np.zeros(shape=(n1,n2,2))
		h2 = np.ones(shape=(n1,n2))
		if polar == True:
			for j in range(n2):
				h2[:,j] = x1
		
		for i in range(n1):
			scrh1 = phi[i,:]
			grad_phi[i,:,1] = self.deriv(scrh1,x2)/h2[i,:]
		for j in range(n2):
			scrh2 = phi[:,j]
			grad_phi[:,j,0] = self.deriv(scrh2,x1)

		return grad_phi

	def Div(self,u1,u2,x1,x2,dx1,dx2,geometry=None):
		(n1, n2) = u1.shape
		Divergence = np.zeros(shape=(n1,n2))
		du1 = np.zeros(shape=(n1,n2))
		du2 = np.zeros(shape=(n1,n2))

		A1 = np.zeros(shape=n1)
		A2 = np.zeros(shape=n2)

		dV1 = np.zeros(shape=(n1,n2))
		dV2 = np.zeros(shape=(n1,n2))

		if geometry == None : geometry = 'cartesian'
		
		#------------------------------------------------
		#  define area and volume elements for the
		#  different coordinate systems
		#------------------------------------------------

		if geometry == 'cartesian' :
			A1[:] = 1.0
			A2[:] = 1.0
			dV1   = np.outer(dx1,A2)
			dV2   = np.outer(A1,dx2)

		if geometry == 'cylindrical' :
			A1 = x1
			A2[:] = 1.0
			dV1 = np.meshgrid(x1*dx1,A2)[0].T*np.meshgrid(x1*dx1,A2)[1].T
			for i in range(n1) : dV2[i,:] = dx2[:]
		
		if geometry == 'polar' :
			A1    = x1
			A2[:] = 1.0
			dV1   = np.meshgrid(x1,A2)[0].T*np.meshgrid(x1,A2)[1].T
			dV2   = np.meshgrid(x1,dx2)[0].T*np.meshgrid(x1,dx2)[1].T

		if geometry == 'spherical' :
			A1 = x1*x1
			A2 = np.sin(x2)
			for j in range(n2): dV1[:,j] = A1*dx1
			dV2   = np.meshgrid(x1,np.sin(x2)*dx2)[0].T*np.meshgrid(x1,np.sin(x2)*dx2)[1].T

		# ------------------------------------------------
		#              Make divergence
		# ------------------------------------------------
		
		
		for i in range(1,n1-1):
			du1[i,:] = 0.5*(A1[i+1]*u1[i+1,:] - A1[i-1]*u1[i-1,:])/dV1[i,:]
		for j in range(1,n2-1):
			du2[:,j] = 0.5*(A2[j+1]*u2[:,j+1] - A2[j-1]*u2[:,j-1])/dV2[:,j]

		print du1[1:10,20]
		Divergence = du1 + du2
		return Divergence

	def curl(self):
		return curlB

	def congrid(self, a, newdims, method='linear', centre=False, minusone=False):
	    '''Arbitrary resampling of source array to new dimension sizes.
	    Currently only supports maintaining the same number of dimensions.
	    To use 1-D arrays, first promote them to shape (x,1).

	    Uses the same parameters and creates the same co-ordinate lookup points
	    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
	    routine of the same name.

	    method:
	    neighbour - closest value from original data
	    nearest and linear - uses n x 1-D interpolations using
				 scipy.interpolate.interp1d
	    (see Numerical Recipes for validity of use of n 1-D interpolations)
	    spline - uses ndimage.map_coordinates

	    centre:
	    True - interpolation points are at the centres of the bins
	    False - points are at the front edge of the bin

	    minusone:
	    For example- inarray.shape = (i,j) & new dimensions = (x,y)
	    False - inarray is resampled by factors of (i/x) * (j/y)
	    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
	    This prevents extrapolation one element beyond bounds of input array.
	    '''
	    if not a.dtype in [np.float64, np.float32]:
		a = np.cast[float](a)

	    m1 = np.cast[int](minusone)
	    ofs = np.cast[int](centre) * 0.5
	    old = np.array( a.shape )
	    ndims = len( a.shape )
	    if len( newdims ) != ndims:
		print "[congrid] dimensions error. " \
		      "This routine currently only support " \
		      "rebinning to the same number of dimensions."
		return None
	    newdims = np.asarray( newdims, dtype=float )
	    dimlist = []

	    if method == 'neighbour':
		for i in range( ndims ):
		    base = np.indices(newdims)[i]
		    dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
				    * (base + ofs) - ofs )
		cd = np.array( dimlist ).round().astype(int)
		newa = a[list( cd )]
		return newa

	    elif method in ['nearest','linear']:
		# calculate new dims
		for i in range( ndims ):
		    base = np.arange( newdims[i] )
		    dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
				    * (base + ofs) - ofs )
		# specify old dims
		olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

		# first interpolation - for ndims = any
		mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
		newa = mint( dimlist[-1] )

		trorder = [ndims - 1] + range( ndims - 1 )
		for i in range( ndims - 2, -1, -1 ):
		    newa = newa.transpose( trorder )

		    mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
		    newa = mint( dimlist[i] )

		if ndims > 1:
		    # need one more transpose to return to original dimensions
		    newa = newa.transpose( trorder )

		return newa
	    elif method in ['spline']:
		oslices = [ slice(0,j) for j in old ]
		oldcoords = np.ogrid[oslices]
		nslices = [ slice(0,j) for j in list(newdims) ]
		newcoords = np.mgrid[nslices]

		newcoords_dims = range(n.rank(newcoords))
		#make first index last
		newcoords_dims.append(newcoords_dims.pop(0))
		newcoords_tr = newcoords.transpose(newcoords_dims)
		# makes a view that affects newcoords

		newcoords_tr += ofs

		deltas = (np.asarray(old) - m1) / (newdims - m1)
		newcoords_tr *= deltas

		newcoords_tr -= ofs

		newa = scipy.ndimage.map_coordinates(a, newcoords)
		return newa
	    else:
		print "Congrid error: Unrecognized interpolation type.\n", \
		      "Currently only \'neighbour\', \'nearest\',\'linear\',", \
		      "and \'spline\' are supported."
		return None





##############--IMAGE ANALYSIS CLASS--##########################
#
# This Class has all the routines for the imaging the data
# and plotting various contours and fieldlines on these images.
# CALLED AFTER pyPLUTO.pload object is defined
#
################################################################



class Image(object):
	def pldisplay(self,var,**kwargs):
		fignum = kwargs.get('fignum',1)
		x1 = kwargs.get('x1')
		x2 = kwargs.get('x2')

		if kwargs.get('polar',False) == True:
			var_cart = self.get_polar_plot(var,**kwargs)
			pcolormesh(var_cart)
			
		
		else:
			if var.shape == var.T.shape  :
				var = var
			else :
				var = var.T

			figure(num=fignum, dpi=80, facecolor='w', edgecolor='k')
			pcolormesh(x1,x2,var,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
		
		title(kwargs.get('title',"Title"),size=kwargs.get('size'))
		xlabel(kwargs.get('label1',"Xlabel"),size=kwargs.get('size'))
		ylabel(kwargs.get('label2',"Ylabel"),size=kwargs.get('size'))
		if kwargs.get('cbar',(False,''))[0] == True:
			colorbar(orientation=kwargs.get('cbar')[1])
	

	def multi_disp(self,*args,**kwargs):
		mvar = []
		var_cart_list=[]
		mfig = figure(num=kwargs.get('fignum',1),facecolor='w', edgecolor='k')
		for arg in args:

			if kwargs.get('polar',False) == True:
				var_cart = self.get_polar_plot(arg,**kwargs)
				var_cart_list.append(var_cart)
			
			else :
				if arg.shape == arg.T.shape  :
					mvar.append(arg)
				else :
					mvar.append(arg.T)
		

		xmin = np.min(kwargs.get('x1'))
		xmax = np.max(kwargs.get('x1'))		
		ymin = np.min(kwargs.get('x2'))
		ymax = np.max(kwargs.get('x2'))

		Ncols = kwargs.get('Ncols')
		Nrows = len(args)/Ncols
		mprod = Nrows*Ncols
		dictcbar=kwargs.get('cbar',(False,'','each'))
		for j in range(mprod):
			mfig.add_subplot(Nrows,Ncols,j+1)
			if kwargs.get('polar',False) == True:
				pcolormesh(var_cart_list[j])
				axis([0.0,var_cart_list[j].shape[1],0.0,var_cart_list[j].shape[0]])
			else:
				imshow(mvar[j],extent=[xmin,xmax,ymin,ymax],origin="image",interpolation="spline36")
			
			xlabel(kwargs.get('label1',mprod*['Xlabel'])[j])
			ylabel(kwargs.get('label2',mprod*['Ylabel'])[j])
			title(kwargs.get('title',mprod*['Title'])[j])
			if (dictcbar[0] == True) and (dictcbar[2] =='each'):
				colorbar(orientation=kwargs.get('cbar')[1])
			if dictcbar[0] == True and dictcbar[2]=='last':
					if (j == np.max(range(mprod))):colorbar(orientation=kwargs.get('cbar')[1])
				
	def field_interp(self,var1,var2,x,y,dx,dy,xp,yp):
		q=[]
		U = var1
		V = var2
		i0 = np.abs(xp-x).argmin()
		j0 = np.abs(yp-y).argmin()
		scrhUx = np.interp(xp,x,U[:,j0])
		scrhUy = np.interp(yp,y,U[i0,:])
		q.append(scrhUx + scrhUy - U[i0,j0])
		scrhVx = np.interp(xp,x,V[:,j0])
		scrhVy = np.interp(yp,y,V[i0,:])
		q.append(scrhVx + scrhVy - V[i0,j0])
		return q

	def field_line(self,var1,var2,x,y,dx,dy,x0,y0):
		xbeg = x[0] - 0.5*dx[0]
		xend = x[-1] + 0.5*dx[-1]

		ybeg = y[0]  - 0.5*dy[0]
		yend = y[-1] + 0.5*dy[-1]

		inside_domain = x0 > xbeg and x0 < xend and y0 > ybeg and y0 < yend


	    ## if inside_domain==False:
	    ##     print " Point %f , %f outside grid range" % (x0,y0)
	    ##     print xbeg, xend
	    ##     print ybeg, yend
	    ##     print np.min(x),np.max(x)
	    ##     print np.min(y),np.max(y)


		MAX_STEPS = 5000
		xln_fwd = [x0]
		yln_fwd = [y0]
		xln_bck = [x0]
		yln_bck = [y0]
		rhs = []
		k = 0

		while inside_domain == True:
		    R1 = self.field_interp(var1,var2,x,y,dx,dy,xln_fwd[k],yln_fwd[k])
		    dl = 0.5*np.max(np.concatenate((dx,dy)))/(np.sqrt(R1[0]*R1[0] + R1[1]*R1[1] + 1.e-14))
		    xscrh = xln_fwd[k] + 0.5*dl*R1[0]
		    yscrh = yln_fwd[k] + 0.5*dl*R1[1]

		    R2 = self.field_interp(var1,var2,x,y,dx,dy,xln_fwd[k],yln_fwd[k])

		    xln_one = xln_fwd[k] + dl*R2[0]
		    yln_one = yln_fwd[k] + dl*R2[1]

		    xln_fwd.append(xln_one)
		    yln_fwd.append(yln_one)
		    inside_domain = xln_one > xbeg and xln_one < xend and yln_one > ybeg and yln_one < yend
		    inside_domain = inside_domain and (k < MAX_STEPS-3)
		    k = k + 1


		k_fwd = k
		qx = np.asarray(xln_fwd[0:k_fwd])
		qy = np.asarray(yln_fwd[0:k_fwd])
		flines={'qx':qx,'qy':qy}
		return flines


	def myfieldlines(self,Data,x0arr,y0arr,stream=False,**kwargs):
		if len(x0arr) != len(y0arr) : print "Input Arrays should have same size"
		QxList=[]
		QyList=[]
		StreamFunction = []
		levels =[]
		if stream == True:
			X, Y = np.meshgrid(Data.x1,Data.x2.T)
			StreamFunction = X*(Data.A3.T)
			for i in range(len(x0arr)):
				nx = np.abs(X[0,:]-x0arr[i]).argmin()
				ny = np.abs(X[:,0]-y0arr[i]).argmin()
				levels.append(X[ny,nx]*Data.A3.T[ny,nx])
			
			contour(X,Y,StreamFunction,levels,colors=kwargs.get('colors'),linewidths=kwargs.get('lw',1),linestyles=kwargs.get('ls','solid'))
		else:
			for i in range(len(x0arr)):
				QxList.append(self.field_line(Data.b1,Data.b2,Data.x1,Data.x2,Data.dx1,Data.dx1,x0arr[i],y0arr[i]).get('qx'))
				QyList.append(self.field_line(Data.b1,Data.b2,Data.x1,Data.x2,Data.dx1,Data.dx1,x0arr[i],y0arr[i]).get('qy'))
				plot(QxList[i],QyList[i],'k--')
			axis([min(Data.x1),max(Data.x1),min(Data.x2),max(Data.x2)])


	def get_polar_plot(self,var,**kwargs):

		def r2theta(outcoords, inputshape, origin):
			xindex, yindex = outcoords
			x0, y0 = origin
			x = xindex - x0
			y = yindex - y0
			r = np.sqrt(x**2 + y**2)
			theta = np.arctan2(np.pi-y, x)
			theta_index = np.round((theta+ np.pi) * inputshape[1] / (1.015*np.pi))
			return (r,theta_index)

		def r2phi(outcoords, inputshape, origin):
			xindex, yindex = outcoords
			x0, y0 = origin
			x = xindex - x0
			y = yindex - y0
			r = np.sqrt(x**2 + y**2)
			phi = np.arctan2(y, x)
			phi_index = np.round((phi + np.pi) * inputshape[1] / (2.1*np.pi))
			return (r,phi_index)

		x1 = kwargs.get('x1') # This is Radial co-ordinate
		x2 = kwargs.get('x2') # This is Either Theta or Phi
		data_arr2d = var

		if kwargs.get('rtheta',False) == True:
			var_cart = S.ndimage.geometric_transform(data_arr2d, r2theta,order=1,output_shape = (data_arr2d.shape[0]*2, data_arr2d.shape[0]),extra_keywords = {'inputshape':data_arr2d.shape,'origin':(data_arr2d.shape[0],0)})

		if kwargs.get('rphi',False) == True:
			var_cart = S.ndimage.geometric_transform(data_arr2d, r2phi, order=1,output_shape = (data_arr2d.shape[0]*2, data_arr2d.shape[0]*2),extra_keywords = {'inputshape':data_arr2d.shape,'origin':(data_arr2d.shape[0],data_arr2d.shape[0])})

		return var_cart
	    
	    
		

	


			
	



      




		    

	    
    
	    
	    
	    
	    
	    
	    



	
        
	  






    
         
    







    
   
