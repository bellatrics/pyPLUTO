import matplotlib
matplotlib.use('TkAgg')



from numpy import arange, sin, pi,log10,max,min,cos
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import pyPLUTO as pp

from Tkinter import *
import sys
import os


class App:
    def __init__(self,master):
        frame = Frame(master)
        
        frame.grid(ipadx=10,ipady=10)
        
        
        self.wdir = os.getcwd() + '/'
        self.I = pp.Image()
        


        self.lb1=Label(frame, text="Nstep").grid(row=0,column=0)
        self.lb2=Label(frame,text="Labels").grid(row=16,column=0)
        
        self.enstep = Entry(frame,width=8)
        self.enstep.grid(row=0,column=1)
        self.enstep.insert(0, "0")

        
        self.varkeys = self.loaddata().get_varinfo(w_dir=self.wdir)['allvars']
        self.grid_dict= self.loaddata().grid(w_dir=self.wdir)
        

        self.ldatabutton=Button(frame,text="Load data",command=self.loaddata)
        self.ldatabutton.grid(row=0,column=2)

        self.ex1 = Entry(frame,width=5)
        self.ex1.grid(row=2,column=0)
        self.ex1.insert(0, "x1")

        self.ex2 = Entry(frame,width=5)
        self.ex2.grid(row=2,column=1)
        self.ex2.insert(0, "x2")
        
        self.ex3 = Entry(frame,width=5)
        self.ex3.grid(row=2,column=2)
        self.ex3.insert(0, "x3")

               

        # place a graph somewhere here
        self.f = Figure(figsize=(8,8), dpi=100)
        self.a = self.f.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.f, master=root)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0,column=3,columnspan=10,rowspan=10,sticky=E)

        



        
        
        self.v = StringVar()
        self.v.set("rho")
        
        for j in range(len(self.varkeys)):
            self.ldata = Radiobutton(frame,text=self.varkeys[j],variable=self.v,value=self.varkeys[j],command=self.getmyvar)
            self.ldata.grid(row=3+j,column=0)
            
        self.slvar = StringVar()
        self.slvar.set("Choose Slice")
        SliceList = ("Along x1","Along x2","Along x3","Along x1-x2","Along x2-x3","Along x3-x1")
        for j in range(len(SliceList)):
            self.sldata = Radiobutton(frame,text=SliceList[j],variable=self.slvar,value=SliceList[j],command=self.setslice)
            self.sldata.grid(row=3+j,column=1)

        #OptionMenu(frame, self.slvar, *SliceList, command=self.setslice).grid(row=1,column=1)

        self.logvar = IntVar()
        self.chkb = Checkbutton(frame,text="Enable Log",variable=self.logvar,onvalue=1,offvalue=0,command=self.logchkcall)
        self.chkb.grid(row=15,column=0)

        self.polarvar = IntVar()
        self.polchkb = Checkbutton(frame,text="Polar",variable=self.polarvar,onvalue=1,offvalue=0,command=self.polchkcall)
        self.polchkb.grid(row=15,column=1)

        self.preaspect = IntVar()
        self.aspectb = Checkbutton(frame,text="Preserve Aspect",variable=self.preaspect,onvalue=1,offvalue=0,command=self.aspchkcall)
        self.aspectb.grid(row=15,column=2)

        self.xlb = Entry(frame,width=15)
        self.xlb.grid(row=16,column=1)
        self.xlb.insert(0, "xlabel")

        self.ylb = Entry(frame,width=15)
        self.ylb.grid(row=16,column=2)
        self.ylb.insert(16, "ylabel")

        self.pltbutton=Button(frame,text="Plot",command=self.plotfinal)
        self.pltbutton.grid(row=18,column=0)

        self.surfbutton=Button(frame,text="Surface",command=self.plotsurface)
        self.surfbutton.grid(row=18,column=1)

        self.clrbutton=Button(frame,text="Clear",command=self.plotclear)
        self.clrbutton.grid(row=18,column=2)

        self.lbinf0 =  Label(frame,text="Information",font=("Times",12,"bold"))
        self.lbinf0.grid(row=25,column=0,sticky=W,columnspan=3)

        self.lbinf1a = Label(frame,text="Dir :",font=("Times",10,"bold")).grid(row=27,column=0,sticky=W,columnspan=3)
        self.lbinf1 =  Label(frame,text=self.wdir).grid(row=28,column=0,sticky=W,columnspan=3)
        self.lbinf2a = Label(frame,text="Domain :",font=("Times",10,"bold")).grid(row=29,column=0,sticky=W,columnspan=3)
        self.lbinf2 = Label(frame,text="n1 x n2 x n3 =  %d x %d x %d " % (self.grid_dict.get('n1'),self.grid_dict.get('n2'),self.grid_dict.get('n3'))).grid(row=30,column=0,sticky=W,columnspan=3)
        self.lbinf3a = Label(frame,text="Time Status",font=("Times",10,"bold")).grid(row=31,column=0,sticky=W,columnspan=3)
        self.lbinf3 = Label(frame,text="Nlast = %d, Time = %f" % (pp.time_info(self.wdir).get('nlast'), pp.time_info(self.wdir).get('time'))).grid(row=32,column=0,sticky=W,columnspan=3)
                        

        
    def loaddata(self):
        mynstep=int(self.enstep.get())
        self.D = pp.pload(mynstep,w_dir=self.wdir)
        return self.D

    def getmyvar(self):
       self.myvar=self.v.get()
        
    def logchkcall(self):
        self.logchk = self.logvar.get()

    def aspchkcall(self):
        self.aspchk=self.preaspect.get()

    def polchkcall(self):
        self.polchk = self.polarvar.get()
        
    def setslice(self):
        self.slicename=self.slvar.get()
        
    def plotclear(self):
        self.a.clear()
       
        
        if len(self.f.axes)>1:
            self.f.delaxes(self.f.axes[1])
            self.f.subplots_adjust(right=0.90)

        self.canvas.show()

    def plotfinal(self):
        if self.logvar.get() == 1:
            self.var = log10(self.D.__getattribute__(self.myvar))
        else:
            self.var = self.D.__getattribute__(self.myvar)

        if self.slicename == "Along x1":
            self.x = self.D.x1
            if self.grid_dict["n3"] == 1:
                self.var = self.var[:,int(self.ex2.get())]
            else:
                self.var = self.var[:,int(self.ex2.get()),int(self.ex3.get())]
            
        elif self.slicename == "Along x2":
            self.x = self.D.x2
            if self.grid_dict["n3"] == 1:
                self.var = self.var[int(self.ex1.get()),:]
            else:
                self.var = self.var[int(self.ex1.get()),:,int(self.ex3.get())]

        else:
            self.x = self.D.x3
            self.var = self.var[int(self.ex1.get()),int(self.ex2.get()),:]

        self.a.set_aspect('auto')
        self.a.plot(self.x,self.var)
        self.a.set_xlabel(self.xlb.get())
        self.a.set_ylabel(self.ylb.get())
        self.canvas.show()

    def plotsurface(self):
       
        if self.logvar.get() == 1:
            self.var = log10(self.D.__getattribute__(self.myvar))
        else:
            self.var = self.D.__getattribute__(self.myvar)
            
        if self.slicename == "Along x1-x2":
            self.x = self.D.x1
            self.y = self.D.x2
            if self.grid_dict["n3"] == 1:
                self.var = self.var[:,:].T
            else:
                if self.polarvar.get() == 1:
                    self.var = self.I.get_polar_plot(self.var[:,:,int(self.ex3.get())],rtheta=True)
                    self.a.axis([0.0,D.n1,0.0,2*D.n1])
                else:
                    self.var = self.var[:,:,int(self.ex3.get())].T
        

        elif self.slicename == "Along x2-x3":
            self.x = self.D.x2
            self.y = self.D.x3
            self.var = self.var[int(self.ex1.get()),:,:].T

        else:
            self.x = self.D.x1
            self.y = self.D.x3
            if self.polarvar.get() == 1:
                self.var = self.I.get_polar_plot(self.var[:,int(self.ex2.get()),:],rphi=True)
                self.a.axis([0.0,2*D.n1,0.0,2*D.n1])
            else:
                self.var = self.var[:,int(self.ex2.get()),:].T

        if self.preaspect.get() == 1:
            self.a.set_aspect('equal')
        else:
            self.a.set_aspect('auto')
            

        if self.polarvar.get() == 1:
           # if self.slicename == "Along x1-x2":self.a.axis([0,max(self.D.n1),0,max(self.D.n2)])
            #if self.slicename == "Along x3-x1":self.a.axis([0,max(self.D.n1),0,max(self.D.n3)])
            self.image=self.a.pcolormesh(self.var,vmin=min(self.var),vmax=max(self.var))
        else:
            self.a.axis([min(self.x),max(self.x),min(self.y),max(self.y)])
            self.image=self.a.pcolormesh(self.x,self.y,self.var,vmin=min(self.var),vmax=max(self.var))
        
        self.a.set_xlabel(self.xlb.get())
        self.a.set_ylabel(self.ylb.get())
        self.f.colorbar(self.image)
        self.canvas.show()
                         

        

    def epssave(self):
        self.f.savefig(self.myvar+'_'+self.enstep.get()+'.eps')
    def pngsave(self):
        self.f.savefig(self.myvar+'_'+self.enstep.get()+'.png')
    def pdfsave(self):
        self.f.savefig(self.myvar+'_'+self.enstep.get()+'.pdf')
    def jpgsave(self):
        self.f.savefig(self.myvar+'_'+self.enstep.get()+'.jpg')
    
    



    
    
            
root=Tk()
app=App(root)
root.title("pyPLUTO")

menubar = Menu(root)
savemenu = Menu(menubar,tearoff=0)
savemenu.add_command(label='EPS',command=app.epssave)
savemenu.add_command(label='PDF',command=app.pdfsave)
savemenu.add_command(label='PNG',command=app.pngsave)
savemenu.add_command(label='JPG',command=app.jpgsave)
menubar.add_cascade(label="Save As", menu=savemenu)



#menubar.add_command(label='Plot',command = app.plotfinal)
#menubar.add_command(label='Surface',command=app.plotsurface)
#menubar.add_command(label='Clear',command=app.plotclear)
menubar.add_command(label='Quit',command=root.quit)

root.config(menu=menubar)

root.mainloop()   


