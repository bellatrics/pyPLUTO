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
        frame.pack()
        
        self.wdir = os.getcwd() + '/'
        self.I = pp.Image()

        
        self.enstep = Entry(frame,width=5)
        self.enstep.pack(side=LEFT)
        self.enstep.insert(0, "0")

        
        self.varkeys = self.loaddata().get_varinfo(w_dir=self.wdir)['allvars']
        self.grid_dict= self.loaddata().grid(w_dir=self.wdir)
        

        self.ldatabutton=Button(frame,text="Load data",command=self.loaddata)
        self.ldatabutton.pack(side=LEFT)
        
        self.v = StringVar()
        self.v.set("rho")
        for variables in self.varkeys:
            self.ldata = Radiobutton(frame,text=variables,variable=self.v,value=variables,command=self.getmyvar)
            self.ldata.pack(side=LEFT)

        self.slvar = StringVar()
        self.slvar.set("Choose Slice")
        SliceList = ("Along x1","Along x2","Along x3","Along x1-x2","Along x2-x3","Along x3-x1")
        OptionMenu(frame, self.slvar, *SliceList, command=self.setslice).pack(side=LEFT)


        self.ex1 = Entry(frame,width=5)
        self.ex1.pack(side=LEFT)
        self.ex1.insert(0, "x1")

        self.ex2 = Entry(frame,width=5)
        self.ex2.pack(side=LEFT)
        self.ex2.insert(0, "x2")
        
        self.ex3 = Entry(frame,width=5)
        self.ex3.pack(side=LEFT)
        self.ex3.insert(0, "x3")

        

       

        self.logvar = IntVar()
        self.chkb = Checkbutton(frame,text="Enable Log",variable=self.logvar,onvalue=1,offvalue=0,command=self.logchkcall)
        self.chkb.pack(side=LEFT)

        self.polarvar = IntVar()
        self.polchkb = Checkbutton(frame,text="Polar",variable=self.polarvar,onvalue=1,offvalue=0,command=self.polchkcall)
        self.polchkb.pack(side=LEFT)

       # self.getplot= Button(frame,text='Plot',fg="black",command=self.plotfinal)
       # self.getplot.pack(side=LEFT)

        # place a graph somewhere here
        self.f = Figure(figsize=(8,8), dpi=100)
        self.a = self.f.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.f, master=root)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        #self.toolbar = NavigationToolbar2TkAgg(self.canvas,frame)
        #self.toolbar.update()
        #self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)


    def loaddata(self):
        mynstep=int(self.enstep.get())
        self.D = pp.pload(mynstep,w_dir=self.wdir)
        return self.D

    def getmyvar(self):
       self.myvar=self.v.get()
        
    def logchkcall(self):
        self.logchk = self.logvar.get()

    def polchkcall(self):
        self.polchk = self.polarvar.get()
        
    def setslice(self,event):
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
                self.var = self.var[int(self.ex2.get()),:,int(self.ex3.get())]

        else:
            self.x = self.D.x3
            self.var = self.var[int(self.ex2.get()),int(self.ex3.get()),:]

        
        self.a.plot(self.x,self.var)
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
            else:
                self.var = self.var[:,int(self.ex2.get()),:].T

        
        
        if self.polarvar.get() == 1:
           # if self.slicename == "Along x1-x2":self.a.axis([0,max(self.D.n1),0,max(self.D.n2)])
            #if self.slicename == "Along x3-x1":self.a.axis([0,max(self.D.n1),0,max(self.D.n3)])
            self.image=self.a.pcolormesh(self.var,vmin=min(self.var),vmax=max(self.var))
        else:
            self.a.axis([min(self.x),max(self.x),min(self.y),max(self.y)])
            self.image=self.a.pcolormesh(self.x,self.y,self.var,vmin=min(self.var),vmax=max(self.var))
        
        
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

menubar = Menu(root)
savemenu = Menu(menubar,tearoff=0)
savemenu.add_command(label='EPS',command=app.epssave)
savemenu.add_command(label='PDF',command=app.pdfsave)
savemenu.add_command(label='PNG',command=app.pngsave)
savemenu.add_command(label='JPG',command=app.jpgsave)
menubar.add_cascade(label="Save As", menu=savemenu)



menubar.add_command(label='Plot',command = app.plotfinal)
menubar.add_command(label='Surface',command=app.plotsurface)
menubar.add_command(label='Clear',command=app.plotclear)
menubar.add_command(label='Quit',command=root.quit)

root.config(menu=menubar)

root.mainloop()   
