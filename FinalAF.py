import numpy as np
from numpy import vectorize
import pickle
#import dill
import time
import types
import random
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import gridspec
from itertools import compress
import math
import copy
def decision(probability):
    return random.random() < probability

gridintime=[] #List for storage of the grid

class heart:
    def __init__(self,L=200,p_unexcitable=0.05,p_fibrosis= 0.8,p_dysf=0.05):

        """#########################PLAN######################
        
      

        1) Don't use a list of refractory cells, therefore they just update the grid with excited cells and then subtract from the array

        2) Use vectors for the excited cells
    

	"""
        self.L=L #Is this the lattice extent (eg length of grid on one axis?) If so, then the grid or edge indices aren't L*L
        self.p_dysf=p_dysf #The fraction of a dysfunctional cells 
        self.p_fibrosis=p_fibrosis #The fraction of missing transversal connections
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
        #self.dysfgrid = np.zeros((self.L, self.L)) 
        #self.edgegrid = np.zeros((self.L, self.L)) + 1
        
        
        """
        n_dysf = (int(self.p_dysf*(self.L**2)))
        cellsdysf=np.random.choice(range(self.L**2),n_dysf,False)
        for elements in cellsdysf:
            self.dysfgrid[int(float(elements)/self.L),elements%self.L]=1
        """

        self.dysfgrid = np.random.rand(self.L, self.L)    #grid of random numbers 
        self.dysfgrid[self.dysfgrid < self.p_dysf] = 1
        self.dysfgrid[self.dysfgrid != 1] = 0
		
	
        self.edgegrid = np.random.rand(self.L, self.L) #grid of random numbers 
        self.edgegrid[self.edgegrid > self.p_fibrosis] = 1
        self.edgegrid[self.edgegrid != 1] = 0
        """
        n_fibrosis = (int(self.p_fibrosis*(self.L**2)))
        fibr=np.random.choice(range(self.L**2),n_fibrosis,False)
        for elements in fibr:
            self.edgegrid[int(float(elements)/self.L),elements%self.L]=0
        """
       
      
    def reinitialise(self):
        
        self.time=0
        self.tcounter=[0]
        self.grid = np.zeros((self.L,self.L))
        self.gridofexcite = copy.deepcopy(self.grid)
        self.dysfgrid = np.random.rand(self.L, self.L) #grid of random numbers 
        self.dysfgrid[self.dysfgrid < self.p_dysf] = 1
        self.dysfgrid[self.dysfgrid != 1] = 0
        
    
        self.edgegrid = np.random.rand(self.L, self.L) #grid of random numbers 
        self.edgegrid[self.edgegrid > self.p_fibrosis] = 1
        self.edgegrid[self.edgegrid != 1] = 0
        """
        self.dysfgrid=self.dysfgrid.flatten()
        np.random.shuffle(self.dysfgrid)
        self.dysfgrid=self.dysfgrid.reshape((self.L,self.L))
        self.edgegrid=self.edgegrid.flatten()
        np.random.shuffle(self.edgegrid)
        self.edgegrid=self.edgegrid.reshape((self.L,self.L))
        """
    
    
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
        
        exciteup = np.roll(self.gridofexcite, 1, axis=0)
        exciteup = (self.grid + exciteup)         
        exciteup[exciteup != self.excitation] = 0   
        exciteup = exciteup*np.roll(self.edgegrid, -1, axis = 0) 
        #exciteup[exciteup != self.excitation] = 0 #only gets true edges
        
        return exciteup
        
    def excitedown(self):
        
        excitedown = np.roll(self.gridofexcite, -1, axis=0) #rolling down
        excitedown = (self.grid + excitedown)             
        excitedown[excitedown != self.excitation] = 0 
        excitedown = excitedown*self.edgegrid
    
        return excitedown
        
    def dysfcheck(self):
        
        randgrid = np.random.rand(self.L, self.L) #grid of random numbers 
        randgrid[randgrid< self.p_unexcitable] = 1
        randgrid[randgrid != 1] = 0
        randgrid = self.dysfgrid*randgrid #gives random number where dysfunctional
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
            self.figure=plt.figure()
            gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1]) 
            self.ax1 = plt.subplot(gs[0], axisbg='black')
            
                
            self.ax1.set_title('Lattice')
            self.ax1.set_xlabel('x')
            self.ax1.set_ylabel('y')
                
            self.ax2 = plt.subplot(gs[1],axisbg='black' )
            self.ax2.set_xlim(0,200)
            self.ax2.set_ylim(-50,50)
               
            
            
            
            self.heart.excitecolumn()
            self.im = self.ax1.imshow(self.heart.grid, extent=(-0,heart.L, heart.L, 0), aspect = 'auto', cmap = "Greys_r")
            self.im.set_cmap('Greys_r') 
            self.im2, = self.ax2.plot(self.heart.tcounter,self.electrocardiot, color='white' )
            
            
            self.interval=1 
            self.counter=0
            if self.replot==False:
                self.anim1 = animation.FuncAnimation(self.figure, self.updatefig,
                            frames=10000, interval=self.interval, blit=True)
            
            if self.store==True:
                #gridintime.append(self.heart.grid) 
                print("self store enacted!")
            
            if self.replot==True:
                print (" elements stored in heart:")
                print (len(gridintime) )
                
                
                self.anim1 = animation.FuncAnimation(self.figure, self.replotfigure,
                            frames=10000, interval=self.interval, blit=True) 
                
        self.timeb=time.time()
        
        print ("Timing",(self.timeb-self.timea))


    """def __getinitargs__(self):
                    temptup2 = (self.heart, self.plot, self.store, self.replot)
                    return temptup2
            
                    #return(self.heart.L, self.heart.p_dysf, self.heart.p_fibrosis, self.heart.p_unexcitable, self.heart.excitation, self.heart.heartbeatsteps, self.heart.grid, self.heart.edges, self.heart.grid_indices, self.heart.edges_indices, self.heart.cells_excited, self.heart.cells_refractory, self.heart.time, self.heart.tcounter, self.heart.cells_dysfindex, self.heart.cells_dysf, self.heart.edges_missingindex, self.heart.edges_missing, self.heart.gridintime)
                
                def __getstate__(self):
                    odict = self.__dict__.copy()
                    return odict  """
    
    def updatefig(self, *args): #Function that yields the data for animation
   
        if (self.heart.time % self.heart.heartbeatsteps)==0 and self.heart.time!=0:    
            self.heart.excitecolumn()

        if self.plot==True and self.store==False:
            self.heart.onestep()
            
        if self.plot==True and self.store==True:
            #self.gridintime.append(self.heart.grid)
            self.heart.onestep()
            
        self.electrocardiot.append(self.heart.electrocardio())
        #self.num_excited.append(len(self.heart.grid[self.heart.grid == self.heart.excitation]))
        
        self.im2.set_data(self.heart.tcounter,self.electrocardiot)
        if len(self.heart.tcounter)>200:
            self.im2.axes.set_xlim(self.heart.tcounter[-200],self.heart.tcounter[-1])
        
        self.im.set_array(self.heart.grid)
        
        self.counter+=1
    
        return [self.im,self.im2]

    def replotfigure(self,*args): #Function to replot stored data
        
        self.im.set_array(self.gridintime[self.counter])   #Updating the grid to animate
        self.counter+=1
        return self.im,   
    
    def storesteps(self):
        
        
        for elements in range(self.stepsstored):
            if (self.heart.time % self.heart.heartbeatsteps)==0 and self.heart.time!=0:    
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

        
        """zero_crossings = np.where(np.diff(np.sign(systems.electrocardiot)))[0][::2] #finds when the electrogram crosses zero and slices to only find odd elements 
                                diff_crossings = np.absolute(zero_crossings-np.roll(zero_crossings, 1)) #to find the time period between crossings that are from +ve to -ve if < 200 then fibrillation starts
                                if any(diff_crossings) < 200:
                                    timefib = zero_crossings[diff_crossings<200]
                                    timefib = len(timefib) 
                                    print("timefib")
                                    print(timefib)"""
         
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
    

h = heart(L=200,p_unexcitable=0.05,p_fibrosis= 0.87,p_dysf=0.05)
#h.electrocardiosetup([100,100])
r = run(heart=h, plot=True,store=False,stepsstored=10000,replot=False)
r.plotecg()
plt.show()"""