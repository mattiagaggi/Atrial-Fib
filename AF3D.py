import numpy as np
from numpy import vectorize
import pickle
#import dill
import time
import types
import random
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import gridspec
import mpl_toolkits
#from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import axes3d, Axes3D 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from itertools import compress
import math
import copy
def decision(probability):
    return random.random() < probability

gridintime=[] #List for storage of the grid

class heart:
    def __init__(self,L=200,p_unexcitable=0.05,p_fibrosis= 0.2,p_dysf=0.05):

        """#########################PLAN######################
        
      

        1) Don't use a list of refractory cells, therefore they just update the grid with excited cells and then subtract from the array

        2) Use vectors for the excited cells
    

	"""
        self.L=L 
        self.p_dysf=p_dysf #The fraction of a dysfunctional cells 
        self.p_fibrosis=p_fibrosis #The fraction of missing transversal connections ( Nu)
        self.p_unexcitable=p_unexcitable #Probablility of a dysfunctional cell being unexcitable
        self.excitation=50 #Value of excitation on the lattice
        self.heartbeatsteps=220 #Time period between excitation wavefronts
        self.grid = np.zeros((self.L,self.L))
        self.gridofexcite = copy.deepcopy(self.grid)
        
        
        self.electrocardiovertpos=np.zeros((self.L,self.L))
        self.electrocardiohorizpos=np.zeros((self.L,self.L))
        self.volt=0
                                        #list of edges True or False for every element self.grid[a,b] the element self.edges[a][b] is the edge below the elements

        self.time=0
        self.tcounter=[0]
   


        self.dysfgrid = np.random.rand(self.L, self.L)    #grid of random numbers 
        self.dysfgrid[self.dysfgrid < self.p_dysf] = 1
        self.dysfgrid[self.dysfgrid != 1] = 0
		
	
        self.edgegrid = np.random.rand(self.L, self.L) #grid of random numbers 
        self.edgegrid[self.edgegrid < self.p_fibrosis] = 1
        self.edgegrid[self.edgegrid != 1] = 0

      
    def reinitialise(self):
        
        self.time=0
        self.tcounter=[0]
        self.grid = np.zeros((self.L,self.L))
        self.gridofexcite = copy.deepcopy(self.grid)

        self.dysfgrid = np.random.rand(self.L, self.L) #grid of random numbers 
        self.dysfgrid[self.dysfgrid < self.p_dysf] = 1
        self.dysfgrid[self.dysfgrid != 1] = 0
         
        self.edgegrid = np.random.rand(self.L, self.L) #grid of random numbers 
        self.edgegrid[self.edgegrid < self.p_fibrosis] = 1
        self.edgegrid[self.edgegrid != 1] = 0

    
    
    def excite(self,a,b): #excitation wavefront
        
        
        self.gridofexcite[a,b] = self.excitation
        self.grid[a,b] = self.excitation
        
        
        
    
    def excitecolumn(self):
        self.gridofexcite[:,0] = self.excitation
        self.grid[:,0] = self.excitation

        
    def onestep(self): #Propagating to the next time step
        
        self.time+=1
        self.tcounter.append(self.time)
        self.repolarisation()
        self.gridofexcite = self.exciteright() + self.exciteleft() + self.exciteup() + self.excitedown() #replaces the list of old excited cells  with the list of newly excited cells
        self.gridofexcite[self.gridofexcite > self.excitation] = self.excitation #making repeated cells = self.excitation
        self.gridofexcite = self.gridofexcite-self.dysfcheck() #getting rid of dysfunctional cells that aren't excited
        self.grid = self.grid + self.gridofexcite # updates the new excited cells on the list of refractory cells
        
      
        
    def repolarisation(self):
        """
        repolarises the grid
        """
        self.grid -=1
        self.grid[self.grid<0] = 0
        
        
    def exciteright(self):
        
        exciteright = np.roll(self.gridofexcite, 1, axis=1) #rolling right
        exciteright[:,0] = 0
        exciteright = (self.grid + exciteright)
        exciteright[exciteright != self.excitation] = 0 #removing refractory cells that would have been excited
        
        return exciteright
    
    def exciteleft(self):

        exciteleft = np.roll(self.gridofexcite, -1, axis=1) #rolling left
        exciteleft[:,self.L-1] = 0
        exciteleft = (self.grid + exciteleft) 
        exciteleft[exciteleft != self.excitation] = 0 
        
        return exciteleft   
    
    def exciteup(self):
        
        exciteup = np.roll(self.gridofexcite, -1, axis=0)
        exciteup = (self.grid + exciteup)         
        exciteup[exciteup != self.excitation] = 0   
        exciteup = exciteup*np.roll(self.edgegrid, -1, axis = 0) 
        #exciteup[exciteup != self.excitation] = 0 #only gets true edges
        
        return exciteup
        
    def excitedown(self):
        
        excitedown = np.roll(self.gridofexcite, 1, axis=0) #rolling down
        excitedown = (self.grid + excitedown)             
        excitedown[excitedown != self.excitation] = 0 
        excitedown = excitedown*self.edgegrid
    
        return excitedown
        
    def dysfcheck(self):
        
        randgrid = np.random.rand(self.L, self.L) #grid of random numbers 
        randgrid[randgrid< self.p_unexcitable] = 1 #Cells that wouldn't be excited, if dysfunctional
        randgrid[randgrid != 1] = 0
        randgrid = self.dysfgrid*randgrid #gets cells that wouldn't be excited and are in the set of dysfunctional cells
        randgrid = self.gridofexcite*randgrid #cells that are dysfunctional and not excited

        
        return randgrid
    
    def Volt(self):
        
        self.volt= 20-((110./self.excitation)*(self.excitation-self.grid))
        
    
        
        
    def electrocardiosetup(self,position):
        
        gridoforizpositions=np.zeros((self.L,self.L))
        gridoforizpositions=(gridoforizpositions+1)*(np.array([range(200)]))
        gridoforizpositions=gridoforizpositions-position[0]
        
        gridofvertpositions=np.zeros((self.L,self.L))
        gridofvertpositions=(gridofvertpositions+1)*(np.transpose(np.array([range(200)])))
        gridofvertpositions=gridofvertpositions-position[1]
        
        self.electrocardiohorizpos=gridoforizpositions
        self.electrocardiovertpos=gridofvertpositions
        
        
    def electrocardio(self):
        self.Volt()
        z=3
        voltxminus1= np.roll(self.volt,1,axis=1)
        voltxminus1[:,0]=0
        voltyminus1= np.roll(self.volt,1,axis=0)
        
        return np.sum ( (self.electrocardiohorizpos *( self.volt-voltxminus1)+self.electrocardiovertpos *(self.volt-voltyminus1))/ ((self.electrocardiovertpos**2 +self.electrocardiohorizpos**2+z**2)**(3./2) ) )

class run: #Class to run code
    def __init__(self,heart,plot=False,store=True,stepsstored=10000,replot=False):
        
        self.heart=heart
        self.timea=time.time()
        self.plot=plot
        self.store=store
        self.replot=replot
        self.gridintime=[]
        self.stepsstored=stepsstored
        self.electrocardiot=[]
        self.num_excited = []

        self.electrocardiot.append(0)
        self.num_excited.append(0)
       
        self.tstartfib=200
        self.tstopfib=210
        self.timecheck=-120 # I set this to be negative so it doesn't fucked up the first time grid[100,100] is excited
        
        
        self.infibrillation=False
        self.tfibrillation=[]
        
        self.fibrillationcounter=0

           
        
        
        
        if self.plot==False and self.replot==False and self.store==True:
            self.timec = time.time()
            self.storesteps()
            self.timed = time.time()
            print("Storage time = %s" %(self.timed-self.timec))
            print("Exclusive self store enacted!")
        
        if self.plot or self.replot:
            fig1 = plt.figure()
            fig2 = plt.figure()
            self.ax1 = fig1.add_subplot(111)
            self.ax1.set_title('Lattice, Nu = %s' %(self.heart.p_fibrosis))
            self.ax1.set_xlabel('x')
            self.ax1.set_ylabel('y')
            self.heart.excitecolumn()
            self.im = self.ax1.imshow(self.heart.grid, extent=(-0,heart.L, heart.L, 0), aspect = 'auto', cmap = "Greys_r")
            self.im.set_cmap('Greys_r')    
            
            
            self.ax = fig2.add_subplot(111, projection='3d')
            self.ax.set_title('3DLattice, Nu = %s' %(self.heart.p_fibrosis))
            self.ax.set_xlabel('x')
            self.ax.set_ylabel('y')
            self.ax.set_zlabel('z')

            print('grid', self.heart.grid)

            u = np.linspace(0, 2 * np.pi, 200)
            v = np.linspace(0, np.pi, 200)

            self.x = 10 * np.outer(np.cos(u), np.sin(v))
            self.y = 10 * np.outer(np.sin(u), np.sin(v))
            self.z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))
            #ax.plot_surface(x, y, z,  color='b')
            #self.surf_plot = self.ax.plot_surface(self.x, self.y, self.z, rstride=1, cstride=1,linewidth = 0, antialiased=False)
            print("yting",len(self.y), len(self.heart.grid))
            facecol = np.absolute(self.heart.grid/self.heart.excitation)
            self.surf = self.ax.plot_surface(self.x, self.y, self.z, rstride=1, cstride=1,linewidth = 0, facecolors=cm.gray(facecol), antialiased=False)
            plt.draw()
            """longitude = x / R
            latitude = 2 * atan(exp(y/R)) - pi/200
            map = Basemap(projection='ortho', lat_0 = 50, lon_0 = -100,
              resolution = 'l', area_thresh = 1000.)
            # plot surface
            map.warpimage(image=self.im)
            # draw the edge of the map projection region (the projection limb)
            map.drawmapboundary()
            # draw lat/lon grid lines every 30 degrees.
            map.drawmeridians(np.linspace(0, 360, 200))
            map.drawparallels(np.linspace(-90, 90, 200))    """
            #matplotlib.interactive(True)
            
            self.interval=1 
            self.counter=0
            if self.replot==False:
                #self.anim1 = animation.FuncAnimation(fig2, self.updatefig,
                #            frames=10000, interval=self.interval, blit=True)
                for i in range(120):
                    self.updatefig()
                    
            
            if self.store==True:
                #gridintime.append(self.heart.grid) 
                print("self store enacted!")
            
            if self.replot==True:
                print (" elements stored in heart:")
                print (len(gridintime) )
                
                
                #self.anim1 = animation.FuncAnimation(self.figure, self.replotfigure,
                #            frames=10000, interval=self.interval, blit=False) 
                
        self.timeb=time.time()
        
        print ("Timing",(self.timeb-self.timea))


    def updatefig(self, *args): #Function that yields the data for animation
   
        if (self.heart.time % self.heart.heartbeatsteps)==0 and self.heart.time!=0:    #why self.heart.time != 0??
            self.heart.excitecolumn()

        if self.plot==True and self.store==False:
            self.heart.onestep()
            
        if self.plot==True and self.store==True:
            self.gridintime.append(self.heart.grid)
            self.heart.onestep()
            
        self.electrocardiot.append(self.heart.electrocardio())
        #self.num_excited.append(len(self.heart.grid[self.heart.grid == self.heart.excitation]))
        
        #self.ax.remove()
        self.surf.remove()
        #self.ax.clear()
        self.im.set_array(self.heart.grid)
        facecol = np.absolute(self.heart.grid/self.heart.excitation)
        self.surf = self.ax.plot_surface(self.x, self.y, self.z, rstride=1, cstride=1,linewidth = 0, facecolors=cm.gray(facecol), antialiased=False)
        #'if surf is not None:
        #ax.collections.remove(surf)
        plt.draw()
        time.sleep(0.01)
        #im3d = map.warpimage(image=self.im)
        self.counter+=1
    

        

    def replotfigure(self,*args): #Function to replot stored data
        
        self.im.set_array(self.gridintime[self.counter])   #Updating the grid to animate
        self.counter+=1
        return self.im,   
    
    def storesteps(self):
        
        
        for elements in range(self.stepsstored):
            if (self.heart.time % self.heart.heartbeatsteps)==0 and self.heart.time!=0:     #why self.heart.time != 0??
                self.heart.excitecolumn()
            self.heart.onestep()
            #self.electrocardiot.append(self.heart.electrocardio())
            self.num_excited.append(len(self.heart.grid[self.heart.grid == self.heart.excitation]))
            

            if len(self.heart.grid[self.heart.grid == self.heart.excitation]) > 220:
                self.fibrillation()

            elif len(self.heart.grid[self.heart.grid==self.heart.excitation]) <= 220 and self.infibrillation == True:
                if len(self.num_excited) > 440 and all(i <= 220 for i in self.num_excited[-441:]) :
                    self.stopfibrillation()

    def fibrillation(self):
            
        if self.infibrillation==False:
            self.infibrillation=True
            self.tfibrillation.append([self.heart.time])
    
    def stopfibrillation(self):

        if self.infibrillation==True:
            #self.fibrillationcounter += 1 #this counts if below fib threshold
            print("in fibrillation")

            self.tfibrillation[-1].append(self.heart.time)
            self.infibrillation=False
            self.fibrillationcounter=0
        
    

    def timeinfibrillation(self):
        timeinfibrillation=0
        print("tfibrillation", self.tfibrillation)
        for elements in self.tfibrillation:
            if len(elements)==2:
                timeinfibrillation+=elements[1]-elements[0]
                print(elements ,"len elements = 2")
            elif len(elements)==1:
                timeinfibrillation+=self.heart.time-elements[0]
                print(elements, "len elements = 1")
                
          
        return timeinfibrillation


         
    def plotecg(self):
        fig3=plt.figure()
        ax = fig3.add_subplot(111)
        plt.plot(self.heart.tcounter,self.num_excited)
        #plt.ylim(-60,60)
        #ax.text(0.1, 0.98,'p_fibros %r' %(self.heart.p_fibrosis), ha='center', va='center',transform=ax.transAxes)
        #ax.text(0.1, 0.94,'p_unex %r' %(self.heart.p_unexcitable), ha='center', va='center',transform=ax.transAxes)
        #ax.text(0.1, 0.9,'p_dysf %r' %(self.heart.p_dysf), ha='center', va='center',transform=ax.transAxes)
        plt.xlabel("time steps")
        plt.ylabel("Voltage")


        
    #def timefib
"""      
plt.show()
   
""" 
h = heart(L=200,p_unexcitable=0.05,p_fibrosis= 0.12,p_dysf=0.05)

r = run(heart=h, plot=True,store=False,stepsstored=10000,replot=False)
#r.plotecg()
plt.show()