import matplotlib
matplotlib.use('TkAgg')



from numpy import arange, sin, pi , cos,log10
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

        self.qbutton= Button(frame,text='Quit',fg="black",command=frame.quit)
        self.qbutton.pack(side=LEFT)
       
        self.clsbutton=Button(frame,text="Clear",command=self.plotclear)
        self.clsbutton.pack(side=LEFT)

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
        SliceList = ("Along x1","Along x2","Along x3")
        OptionMenu(frame, self.slvar, *SliceList, command=self.setslice).pack(side=LEFT)

        
        
        self.ex3 = Entry(frame,width=5)
        self.ex3.pack(side=LEFT)
        self.ex3.insert(0, "x3")

        self.ex2 = Entry(frame,width=5)
        self.ex2.pack(side=LEFT)
        self.ex2.insert(0, "x2")

        self.ex1 = Entry(frame,width=5)
        self.ex1.pack(side=LEFT)
        self.ex1.insert(0, "x1")

        self.logvar = IntVar()
        self.chkb = Checkbutton(frame,text="Enable Log",variable=self.logvar,onvalue=1,offvalue=0,command=self.logchkcall)
        self.chkb.pack(side=LEFT)

        self.getplot= Button(frame,text='Plot',fg="black",command=self.plotfinal)
        self.getplot.pack(side=LEFT)

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
        
    def setslice(self,event):
        self.slicename=self.slvar.get()
        
    def plotclear(self):
        self.a.clear()
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
    

    
    
            
root=Tk()
app=App(root)
root.mainloop()   
