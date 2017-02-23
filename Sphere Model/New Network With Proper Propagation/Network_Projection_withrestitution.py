import numpy as np
from scipy.sparse import csr_matrix
import operator
from itertools import chain
import time
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
from matplotlib import collections  as mc
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import types
import pickle
import ico_prop_projection as sp

def decision(p):
    if np.random.random()<p:
        return True
    else:
        return False



def trial():
    nodes=range(40000)
    a=np.arange(40000)
    b=np.arange(40000)
    c=np.arange(40000)
    np.random.shuffle(b)
    np.random.shuffle(c)
    l=[]
    for e in np.arange(40000):
        l.append((a[e],b[e]))
        l.append((c[e],b[e]))
        l.append((a[e],c[e]))
    array_transv=[]
    p=[0,1]
    return nodes,l,array_transv,1,p,0.5,0.5
        
    

class create_network:
    def __init__(self,array_nodesindices,array_vertical,array_transv,p_transv,impulse_start,p_dysf,p_unexcitable, excitation, hbs,recursion=5):
        """
        this inputs array_nodesindices as a list of faces indices,
        
        array_vertical and array_transv as  list of tuples
        """
        t0=time.time()
            
        # watch it because here we assume that faces go from 0 to n
        
        self.p_dysf=p_dysf
        self.p_transv=p_transv
        self.p_unexcitable=p_unexcitable
        self.excitation=excitation
        self.heartbeatssteps=hbs
        self.totalruns = 0
        
        self.recursion=recursion
        print ("recursion=",self.recursion,"make sure the restitution curve is set accordingly")
        
        
        self.array_vertical=array_vertical
        self.impulse_start=impulse_start
        self.size=len(array_nodesindices)
        self.nodes=np.zeros(self.size)
       
        self.array_transv=[]
        self.excited=[]
        self.dysf=np.random.rand(self.size)  #create array of dysf cells
        self.dysf[self.dysf < self.p_dysf] = 1
        self.dysf[self.dysf != 1] = 0
        self.unexcited=np.random.rand(self.size)
        self.array_alltransv=array_transv
       
        
        for elements in self.array_alltransv: #append transversal connections to a new list according to the probability of transv connect
          if decision(self.p_transv):
              self.array_transv.append(elements)
        
        self.connections = {}   #create dictionary with connections ex d={1:[0],2:[3,4],.....}
        for key,value in (self.array_vertical+self.array_transv):
            try:
                self.connections[key].append(value)
            except KeyError:
                self.connections[key]=[value]
        
        for value,key in (self.array_vertical+self.array_transv): #append also connections from second to first elements
            try:
                self.connections[key].append(value)
            except KeyError:
                self.connections[key]=[value]
        
        #print ('initialisation time:',(t1-t0))
     
        self.refr=np.zeros(self.size)
        self.tempgrid=np.zeros(self.size)
        self.tempgrid+=220
        
       
        
        if self.recursion==6:
            a=2
        elif self.recursion==5:
            a=4
        #in the case of the sphere we coarse grain the refractory period by an additional factor depending on the recursion
        #for recursion 5 instead of 50 we have 12.5 therefore the factor a is 4
        #for recursion 6 .......................25 ...........................2
        
        x=np.arange(0,800,3*a) #this is the time between impulses we choose step a*3 because we will divide by a*3 later
        y=np.array([])
        
        #change this
        for elements in x:   #make vector y with actual APD
            if elements<=210:
                y=np.append(y,90)
            if elements>210 and  elements<=400:
                y=np.append(y, ((130./190)*elements+(220-400*130/190)))
            if elements>400:
                y=np.append(y,220)
        
        z=y*150./220    #make z which is steps
        stepsx=x/(3.*a) #coarse grain by a*3
        stepsx=np.round(stepsx) #we need to round
        stepsx=stepsx.astype(int)
        stepsz=z/(3*a)
        stepsz=np.round(stepsz)
        
        while len(stepsz)<40000:
            
            stepsz=np.append(stepsz,(50./a))
       
        stepsz=stepsz.astype(int)      
            
        self.maparray=stepsz
        t1=time.time()
        print("network creation time", t1-t0)
        
    
    
    
    
    def excite(self,nodeindex):
        

     
        self.excited.append(nodeindex)
        self.nodes[nodeindex]=self.excitation

    def onestep(self):

        

        newnodes=np.zeros(self.size)
       
        self.tempgrid+=1
        
        if len(self.excited)!=0:
            f = operator.itemgetter(*self.excited) # works well with list or array
            
            
            self.excited=f(self.connections) # alternative function: list(self.connections[_] for _ in self.excited) check which one is faster
            try:
                self.excited=list(chain(*self.excited))#the output of f() is ( [....],[...],....) this changes into a unique list
            except TypeError:
                self.excited=list(chain(self.excited))
            
        
            newnodes[self.excited]+=1   #this seems to take more than twice as much the previous operation 0.002
            #need to account for repeated elements
            
            self.unexcited=np.random.rand(self.size)  #create array of dysf cells maybe can find a quicker way #this seems to take long 0.001
            self.unexcited[self.unexcited < self.p_unexcitable] = 1
            self.unexcited[self.unexcited != 1] = 0
           
            
            newnodes[newnodes>0]=1
            newnodes=newnodes-(self.dysf*self.unexcited*newnodes) #remove excited dysf cells which  are not excited
            newnodes*= (self.nodes==0) #removes refractory cells
        
            self.excited=np.flatnonzero(newnodes)
     
            self.excited=self.excited.tolist()
    
      
        self.nodes-=1  
        self.nodes[self.nodes==-1]=0
        
        #maybe switch these two
        self.refr=self.maparray[self.tempgrid.astype(int)]
        
      
        
        self.nodes += newnodes*self.refr
        
        self.tempgrid=self.tempgrid*(newnodes==0)
        #print ("excited cells are",self.excited)
     
        
    
   
   
    
    def reinitialise(self):
        self.nodes=np.zeros(self.size)
       
        self.array_transv=[]
        self.excited=[]
        self.dysf=np.random.rand(self.size)  #create array of dysf cells
        self.dysf[self.dysf < self.p_dysf] = 1
        self.dysf[self.dysf != 1] = 0
        self.unexcited=np.random.rand(self.size)
        
        self.refr=np.zeros(self.size)
        self.tempgrid=np.zeros(self.size)
        self.tempgrid+=220
                
        for elements in self.array_alltransv: #append transversal connections to a new list according to the probability of transv connect
          if decision(self.p_transv):
              self.array_transv.append(elements)
         
        self.connections = {}   #create dictionary with connections ex d={1:[0],2:[3,4],.....}
        for key,value in (self.array_vertical+self.array_transv):
            try:
                self.connections[key].append(value)
            except KeyError:
                self.connections[key]=[value]
        
        for value,key in (self.array_vertical+self.array_transv): #append also connections from second to first elements
            try:
                self.connections[key].append(value)
            except KeyError:
                self.connections[key]=[value]
      
        


class run:
    
    def __init__(self,network, plot=False,store=True,runs=5000, fib_threshold=350):
        
        
        self.runs=runs   
        self.plot= plot
        self.store=store
        self.network=network
        self.sph = []
        self.time=0
        self.nodeshistory=[]
        self.excitedhistory=[]
        self.nodeshistory.append(self.network.nodes)
        self.num_excited=[]
        
        self.fib_threshold=fib_threshold
        self.infibrillation=False
        self.tfibrillation=[]
        
        self.zeroexcited=[]
        
        
    def propagate_storage(self):
      
        #print ("you set store==", self.store)
        for times in range(self.runs):
                
            if self.network.totalruns%self.network.heartbeatssteps == 0: #self.time%self.network.heartbeatssteps==0:
                for elements in self.network.impulse_start:
                    self.network.excite(elements)
                    
            self.time+=1
            self.network.totalruns+=1
            if self.store==True:
                #self.excitedhistory.append(self.network.excited)
                self.num_excited.append(len(self.network.excited))
                
                if self.num_excited[-1]>=self.fib_threshold: # if more excited cells than threshold enters fib
                    self.fibrillation()
                elif self.infibrillation==True:   # if in fibrillation  and if the num excited has been below the threshold for two cicles(heartbeatsteps) then it stops fib
                    if len(self.num_excited) > 2*self.network.heartbeatssteps:
                         if all(i <= self.fib_threshold for i in self.num_excited[-(1+2*self.network.heartbeatssteps):]):
                             self.stopfibrillation()
                    
                
                
            self.network.onestep()

    
   
        
    def fibrillation(self):
            
        if self.infibrillation==False:
            self.infibrillation=True
            self.tfibrillation.append([self.time])
    
    
    def stopfibrillation(self):

        self.tfibrillation[-1].append(self.time)
        self.infibrillation=False
        
        
    

    def timeinfibrillation(self):
        timeinfibrillation=0
        print("tfibrillation", self.tfibrillation)
        for elements in self.tfibrillation:
            if len(elements)==2:
                timeinfibrillation+=elements[1]-elements[0]
                print(elements ,"len elements = 2")
            elif len(elements)==1:
                timeinfibrillation+=self.time-elements[0]
                print(elements, "len elements = 1")
                
          
        return timeinfibrillation
    
    def animator(self,sph):
        
        print( "you set store==",self.store)
        self.sph = sph

        self.surf, self.fig = self.plot_sphere_a(sph)
        self.sph = sph
        self.anim1 = animation.FuncAnimation(self.fig, self.updatefig, fargs = (colours, self.surf),
                frames=self.runs, interval=25, blit=False)
        if self.store==True:
                self.excitedhistory.append(self.network.excited)
                self.num_excited.append(len(self.network.excited))
        plt.show()

          
            
                
                   
                
                
    def propagate_a(self):
        for times in range(runs):
            if self.store==True:
                
                if self.network.totalruns%self.network.heartbeatssteps == 0: #self.time%self.network.heartbeatssteps==0:
                    for elements in self.network.impulse_start:
                        self.network.excite(elements)
                    
                self.time+=1
                self.network.totalruns+=1
                self.excitedhistory.append(self.network.excited)
            self.network.onestep()
            colours = self.network.nodes
            colours = colours/sum(colours) 
            yield colours
    
            
                
    def long_edges(x, y, triangles, radio=22):
        out = []
        for points in triangles:
            #print points
            a,b,c = points
            d0 = np.sqrt( (x[a] - x[b]) **2 + (y[a] - y[b])**2 )
            d1 = np.sqrt( (x[b] - x[c]) **2 + (y[b] - y[c])**2 )
            d2 = np.sqrt( (x[c] - x[a]) **2 + (y[c] - y[a])**2 )
            max_edge = max([d0, d1, d2])
            #print points, max_edge
            if max_edge > radio:
                out.append(True)
            else:
                out.append(False)
        return out               
                
                
    def plot_sphere_a(self, sph):
        self.fig = plt.figure()
        #self.ax = self.fig.add_subplot(111, projection='3d')#.gca(projection='3d')projection = 'mollweide'
        #self.ax = self.fig.add_subplot(111)#, projection = '3d') #'mollweide')#.gca(projection='3d')projection = 'mollweide'
        #self.ax.view_init(elev=90., azim=0)
        
        #self.ax.axis([-1,1,-1,1, -1, 1])
        """
        self.ax.set_xlim(-0.55, 0.55)
        self.ax.set_ylim(-0.55, 0.55)
        self.ax.set_zlim(-0.55, 0.55)
        """
        #plt.gca().set_aspect('equal')
        #axes = plt.gca()
        #axes.set_xlim([-4,4])
        #axes.set_ylim([-5,35])
        self.x, self.y, self.z   = sph.ch.points[sph.ch.vertices][:,0],sph.ch.points[sph.ch.vertices][:,1], sph.ch.points[sph.ch.vertices][:,2]
        #print("self.x", self.x)
        #self.avg = (self.xc+ self.yc+ self.zc)/3
        #self.x, self.y = sph.Mercator_Projection(-2.3*np.pi/8.)#sph.Mercator_Projection()#s#Mollewide
        self.z = np.zeros(self.x.shape)
        self.triangles = sph.ch.simplices
        #self.xy=np.vstack((self.x, self.y)).T
        #tri = Delaunay(self.xy)
        #tris =mtri.Triangulation(self.x, self.y)
        #self.triangles = tris.triangles
        #self.ax.tricontour(triangles = self.triangles, Z=self.z)#, cmap=plt.cm.Greys_r, antialiased=False)
        #self.surf = self.ax.
        #self.triangles = tri.simplices
        #mask = [s]
        self.surf =plt.tripcolor(self.x, self.y,  self.triangles, facecolors = self.triangles[:,0])#facecolours = self.triangles[:,0], edgecolors = 'k')#, cmap=plt.cm.Greys_r, antialiased=False)
        #self.surf = self.ax.plot_trisurf(self.x,self.y,self.z, triangles=sph.ch.simplices, cmap=plt.cm.Greys_r, alpha = 1)
        #self.surf.set_array(colours)
        #sph.icosahedron_vertices=np.asarray(sph.icosahedron_vertices)
        #self.ax.scatter(sph.icosahedron_vertices[:,0],sph.icosahedron_vertices[:,1], sph.icosahedron_vertices[:,2], c='red')
        
        plt.show()
        return self.surf, self.fig
        

    def updatefig(self, *args): #Function that yields the data for animation
        
        if self.network.totalruns%self.network.heartbeatssteps == 0: #self.time%self.network.heartbeatssteps==0:
            for elements in self.network.impulse_start:
                self.network.excite(elements)
                
        self.time+=1
        self.network.totalruns+=1
        if self.store==True:
            self.excitedhistory.append(self.network.excited)
        self.network.onestep()
        colours = self.network.nodes
        colours = colours/sum(colours) 
        
        print("length of colours and x",len(colours), len(self.x))
        #
        #self.ax.clear()
        self.surf = plt.tripcolor(self.x, self.y,  self.triangles, self.z, facecolors = colours[:len(self.triangles)], cmap=plt.cm.Greys_r, antialiased=False)
        #self.surf = self.ax.plot_trisurf(self.x,self.y,self.z, triangles=self.triangles, cmap=plt.cm.Greys_r, antialiased=False)
        #self.surf.set_array(colours)
        """
        self.ax.set_xlim(-0.55, 0.55)
        self.ax.set_ylim(-0.55, 0.55)
        self.ax.set_zlim(-0.55, 0.55)
        """
        print("time", self.time)
        return self.surf,
                

class Define_Connections:
    
    
    def __init__(self,s):
        self.s = s
        
    def num_func(self, recursion_lvl):
        t_0 = 3
        for x in range(recursion_lvl-2):
            t_n = 2*t_0 + 3
            t_0 = t_n
        return t_n
       
    def define_connections(self):
        
        self.s.construct_icosphere()
        vertex1 = self.s.ch.points[self.s.ch.simplices[:,0]]
        self.colours = vertex1
        pent_faces, pent_ind = self.s.find_pentagon(ind = 0)
        
        
        next_tri, vconn = self.s.next_row_tri_v(pent_ind, pent_ind) #Defines vertical connections 
        face_cache = np.hstack((pent_ind, next_tri)) 
        next_tri2, hconn, x_v_conn=  self.s.next_row_tri_h(next_tri, face_cache) 
        vconn = vconn + x_v_conn

        self.startimp = np.hstack((face_cache, next_tri2))
           
        rang_val = int((self.num_func(self.s.recursion_level)+1)/2.)  #-2
        for i in range(rang_val):
            face_cache = np.hstack((face_cache, next_tri2))
            next_tri3, vconn2 =  self.s.next_row_tri_v(next_tri2, face_cache)
            face_cache = np.hstack((face_cache, next_tri3))
            next_tri4, hconn2 , x_v_conn=  self.s.next_row_tri_h(next_tri3, face_cache)
            next_tri2 = next_tri4
            vconn = vconn + vconn2 + x_v_conn
            hconn = hconn + hconn2
        
        face_cache = np.hstack((face_cache, next_tri2))
        next_tri3, vconn2 =  self.s.next_row_tri_v(next_tri2, face_cache)
        vconn = vconn + vconn2
        
        
        
        pent_faces2, pent_ind2 = self.s.find_pentagon(ind = 3)
        next_trif, vconnf = self.s.next_row_tri_v(pent_ind2, pent_ind2) #Defines vertical connections 
        face_cachef = np.hstack((pent_ind2, next_trif)) 
        next_tri2f, hconnf, x_v_connf=  self.s.next_row_tri_h(next_trif, face_cachef) 
        vconnf = vconnf + x_v_connf
           
        rang_val = int((self.num_func(self.s.recursion_level)+1)/2.)
        for i in range(rang_val):
            face_cachef = np.hstack((face_cachef, next_tri2f))
            next_tri3f, vconn2f =  self.s.next_row_tri_v(next_tri2f, face_cachef)
            face_cachef = np.hstack((face_cachef, next_tri3f))
            next_tri4f, hconn2f , x_v_connf=  self.s.next_row_tri_h(next_tri3f, face_cachef)
            next_tri2f = next_tri4f
            vconnf = vconnf + vconn2f + x_v_connf
            hconnf = hconnf + hconn2f
        
        face_cachef = np.hstack((face_cachef, next_tri2f))
        next_tri3f, vconn2f =  self.s.next_row_tri_v(next_tri2f, face_cachef)
        vconnf = vconnf + vconn2f
        
        vconn = vconn + vconnf
        hconn = hconn + hconnf
        self.vconn = self.s.remove_duplicates(vconn)
        self.hconn = self.s.remove_duplicates(hconn)
        
        return self.colours, self.vconn,self.hconn, self.startimp#pent_ind
    
    
    def connectome(self):
        vertex1, vertex2, vertex3 = self.s.ch.points[self.s.ch.simplices[:,0]],self.s.ch.points[self.s.ch.simplices[:,1]], self.s.ch.points[self.s.ch.simplices[:,2]]
        phi_avg = (vertex1+vertex2+vertex3)/3.
        
        linesv = [(tuple(phi_avg[i[0]]), tuple(phi_avg[i[1]]) )for i in self.vconn]
        colv = [(1,0,0,1) for i in range(len(linesv))]
        linesh = [(tuple(phi_avg[i[0]]), tuple(phi_avg[i[1]])) for i in self.hconn]
        colh = [(0,0,1,1) for i in range(len(linesh))]
        colx = [(0,1,0,1) for i in range(len(linesh))]
        coly = [(1,1,0,1) for i in range(len(linesh))]
        
        lines = linesv+linesh
        c = np.array(colv+colh)
    
        
        
        lc1 = Line3DCollection(linesv[int(len(linesv)/2):], colors=colv, linewidths=3)#[:int(5*len(linesv)/7)] [:int(len(linesv)/2)] + linesv[int(3*len(linesv)/4):]
        lc2 = Line3DCollection(linesh[int(len(linesh)/2):], colors=colh, linewidths=3)#[:int(5*len(linesh)/7)]
        #lc3 = Line3DCollection(linesv[int(len(linesv)/2):], colors=colx, linewidths=3)
        #lc4 = Line3DCollection(linesh[int(len(linesh)/2):], colors=coly, linewidths=3)
        #fig, ax = pl.subplots()
        #ax.add_collection(lc)
        fig2 = plt.figure()
        ax2 = fig2.gca(projection='3d')
        #ax2.plot(phi_avg[:,0], phi_avg[:,1], phi_avg[:,2], 'r', markersize=1)
        ax2.add_collection(lc1)
        ax2.add_collection(lc2)
        plt.show()



"""
#opening pickled files for reinstatement of connections
s.construct_icosphere()

d=open('horiz_conn_rec_6xconn.pkl', 'rb')
hconn=pickle.load(d)

f=open('vert_conn_rec_6xconn.pkl', 'rb')
vconn=pickle.load(f)

e=open('startimp_ind_rec_6xconn.pkl', 'rb')
pent_ind=pickle.load(e)

g=open('colours_rec_6xconn.pkl', 'rb')
colours=pickle.load(g)

#h=open('sph_rec_6xconn.pkl', 'rb')
#s=pickle.load(h)

"""

"""
s = sp.Sphere( recursion_level = 6 )
conn = Define_Connections(s)
colours, vconn, hconn, pent_ind = conn.define_connections() #not needed if using pickled data



n = create_network(array_nodesindices = np.arange(len(colours)),
                   array_vertical = vconn,
                   array_transv = hconn,
                   p_transv = 0.1,
                   impulse_start = pent_ind,
                   p_dysf = 0.1,
                   p_unexcitable = 0.1,
                   excitation = 12, 
                   hbs = 55)

runc = run(network = n, plot=True,store=True,runs=1000)

"""
#for storing data instead

#runc = run(network = n, plot=False,store=True,runs=1000,fib_threshold=322)
#runc.propagate_storage()
#runc.animator(s)

"""
for recursion 5 
the path from the top to the bottom is 47
and the width of the equator 150 

instead of coarse graining the model by 5 as in Kims we coarse grain it by 20 (we have 47 instead of 200)
this gives refractory period equal to 13(12.5) and excitation of 55( using T=660ms)
this gives refractory period equal to 13 and excitation of 63( using T=760ms)

180 fib threshold

"""

"""
for recursion 6 
the path from the top to the bottom is 97
and the width of the equator 320

instead of coarse graining the model by 5 as in Kims we coarse grain it by 10 (we have 97 instead of 200)
this gives refractory period equal to 25 and excitation of 110( using T=660ms)
this gives refractory period equal to 25 and excitation of 126( using T=760ms)

fib threshold 350

"""




#conn.connectome() #visualisation of connections, as you've seen


#############################
#Storing animation 

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=25, metadata=dict(artist='Me'), bitrate=1800)
#runc.anim1.save('Recursion6xvom.mp4', fps=25, extra_args=['-vcodec', 'libx264'])
#runc.anim1.save('')
#plt.show()
#############################


"""
#closing pickle files
d.close()
e.close()
f.close()
g.close()
#h.close()
"""

"""
ONLY NEEDED WHEN PICKLING NEW DATA
fileobj1 = open('horiz_conn_rec_6xconn.pkl', 'wb')
pickle.dump(hconn, fileobj1, -1)

fileobj2 = open('vert_conn_rec_6xconn.pkl', 'wb')
pickle.dump(vconn, fileobj2, -1)

fileobj0 = open('startimp_ind_rec_6xconn.pkl', 'wb')
pickle.dump(pent_ind, fileobj0, -1)

fileobj3 = open('colours_rec_6xconn.pkl', 'wb')
pickle.dump(colours, fileobj3, -1)

fileobjsph = open('sph_rec_6xconn.pkl', 'wb')
pickle.dump(s, fileobjsph, -1)


fileobj1.close()
fileobj2.close()
fileobj0.close()
fileobj3.close()
fileobjsph.close()
"""


