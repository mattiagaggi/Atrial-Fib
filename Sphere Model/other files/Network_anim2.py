import numpy as np
import operator
from itertools import chain
import time
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import ico_propagation as sp


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
    def __init__(self,array_nodesindices,array_vertical,array_transv,p_transv,impulse_start,p_dysf,p_unexcitable):
        """
        this inputs array_nodesindices as a list of faces indices,
        
        array_vertical and array_transv as  list of tuples
        """
        t0=time.time()
            
        # watch it because here we assume that faces go from 0 to n
        
        self.p_dysf=p_dysf
        self.p_transv=p_transv
        self.p_unexcitable=p_unexcitable
        self.excitation=50
        self.heartbeatssteps=220
        self.totalruns = 0
    
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
        t1=time.time()
        print ('initialisation time:',(t1-t0))
     
    def excite(self,nodeindex):
        

     
        self.excited.append(nodeindex)
        self.nodes[nodeindex]=self.excitation

    def onestep(self):

        

        newnodes=np.zeros(self.size)
        time1=time.time()
    
        
        
        if len(self.excited)!=0:
            f = operator.itemgetter(*self.excited) # works well with list or array
            time2=time.time()
            print("time taken",(time2-time1))
            
            
            self.excited=f(self.connections) # alternative function: list(self.connections[_] for _ in self.excited) check which one is faster
            try:
                self.excited=list(chain(*self.excited))#the output of f() is ( [....],[...],....) this changes into a unique list
            except TypeError:
                self.excited=list(chain(self.excited))
            time3=time.time()
            print("time taken",(time3-time2)) #seems quick but can improve
        
            newnodes[self.excited]+=1   #this seems to take more than twice as much the previous operation 0.002
            #need to account for repeated elements
            time4=time.time()
            print("time taken",(time4-time3))
        
            self.unexcited=np.random.rand(self.size)  #create array of dysf cells maybe can find a quicker way #this seems to take long 0.001
            self.unexcited[self.unexcited < self.p_unexcitable] = 1
            self.unexcited[self.unexcited != 1] = 0
            time5=time.time()
            print("time taken",(time5-time4))
            
            newnodes[newnodes>0]=1
            newnodes=newnodes-(self.dysf*self.unexcited*newnodes) #remove excited dysf cells which  are not excited
            newnodes*= (self.nodes==0) #removes refractory cells
        
            self.excited=np.flatnonzero(newnodes)
            self.excited=self.excited.tolist()
    
        time6=time.time()
        self.nodes-=1  
        self.nodes[self.nodes==-1]=0
        self.nodes += newnodes*self.excitation
        print (np.flatnonzero(self.nodes))
        time7=time.time()
        print("time taken",(time7-time6))
        
        print ("total time", ( time7-time1))
   
    def reinitialise(self):
        self.nodes=np.zeros(self.size)
       
        self.array_transv=[]
        self.excited=[]
        self.dysf=np.random.rand(self.size)  #create array of dysf cells
        self.dysf[self.dysf < self.p_dysf] = 1
        self.dysf[self.dysf != 1] = 0
        self.unexcited=np.random.rand(self.size)
                
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
    
    def __init__(self,network, plot=False,store=True,runs=5000):
        
        
        self.runs=runs   
        self.plot= plot
        self.store=store
        self.network=network
        self.sph = []
        self.time=0
        self.nodeshistory=[]
        self.excitedhistory=[]
        self.nodeshistory.append(self.network.nodes)
        
        
    def propagate_n(self):
        for times in range(self.runs):
            if self.store==True:
                
                if self.network.totalruns%self.network.heartbeatssteps == 0: #self.time%self.network.heartbeatssteps==0:
                    for elements in self.network.impulse_start:
                        self.network.excite(elements)
                    
                self.time+=1
                self.network.totalruns+=1
                self.excitedhistory.append(self.network.excited)
            self.network.onestep()
            
    def propagate_a(self):
        for times in range(self.runs):
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
    
    def animator(self,sph):
        self.sph = sph
        if self.plot == True:
            self.surf, self.fig = self.plot_sphere_a(sph)
            self.sph = sph
            colours = self.propagate_a()
            #colours = self.network.nodes
            #colours = colours/sum(colours) 
            self.anim1 = animation.FuncAnimation(self.fig, self.updatefig, fargs = (colours, self.surf),
                    frames=self.runs, interval=1, blit=False)
            #for i in range(self.runs):
            #    self.updatefig(sph)
            plt.show()
          
            
                
                   
                
                
    def plot_sphere_a(self, sph):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')#.gca(projection='3d')
        self.x, self.y, self.z   = sph.ch.points[sph.ch.vertices][:,0],sph.ch.points[sph.ch.vertices][:,1], sph.ch.points[sph.ch.vertices][:,2]
        self.surf = self.ax.plot_trisurf(self.x,self.y,self.z, triangles=sph.ch.simplices, cmap=plt.cm.Greys_r)
        self.triangles = sph.ch.simplices
        sph.icosahedron_vertices=np.asarray(sph.icosahedron_vertices)
        self.ax.scatter(sph.icosahedron_vertices[:,0],sph.icosahedron_vertices[:,1], sph.icosahedron_vertices[:,2], c='red')
        
        plt.show()
        return self.surf, self.fig
    
        """    
        if self.plot or self.replot:
            fig1 = plt.figure()
            fig2 = plt.figure() 
            
            
            self.ax = fig2.add_subplot(111, projection='3d')
            self.ax.set_title('3DLattice, Nu = %s' %(1-self.network.p_transv))
            self.ax.set_xlabel('x')
            self.ax.set_ylabel('y')
            self.ax.set_zlabel('z')

           
            facecol = colours
            self.surf = self.ax.plot_surface(self.x, self.y, self.z, rstride=1, cstride=1,linewidth = 0, facecolors=cm.gray(facecol), antialiased=False)
            #plt.draw()
            
            self.interval=25
            self.counter=0
        """
        

    def updatefig(self, *args): #Function that yields the data for animation
        
        if self.network.totalruns%self.network.heartbeatssteps == 0: #self.time%self.network.heartbeatssteps==0:
            for elements in self.network.impulse_start:
                self.network.excite(elements)
                
        self.time+=1
        self.network.totalruns+=1
        self.excitedhistory.append(self.network.excited)
        self.network.onestep()
        colours = self.network.nodes
        colours = colours/sum(colours) 
        
        
        self.ax.clear()
        self.surf = self.ax.plot_trisurf(self.x,self.y,self.z, triangles=self.triangles, cmap=plt.cm.Greys_r, antialiased=False)
        self.surf.set_array(colours)
        
        return self.surf,
        
                




#def conn_4():
#    s = sp.Sphere(vertices = icosahedron_vertices, faces = faces, recursion_level = 4 )
#    
#    s.construct_icosphere()
#    x, y, z   = s.ch.points[s.ch.vertices][:,0],s.ch.points[s.ch.vertices][:,1], s.ch.points[s.ch.vertices][:,2]
#    
#    
#    vertex1, vertex2, vertex3 = s.ch.points[s.ch.simplices[:,0]],s.ch.points[s.ch.simplices[:,1]], s.ch.points[s.ch.simplices[:,2]]# hull_pts[ch.simplices[:,1]], hull_pts[ch.simplices[:,2]]
#    
#    
#    
#    #CHANGE TO SPHERICAL POLARS
#    sph_v_1 = s.sph_polar_convert(vertex1)
#    sph_v_2 = s.sph_polar_convert(vertex2) 
#    sph_v_3 = s.sph_polar_convert(vertex3) 
#    
#    phi_avg = (sph_v_1[:,2] + sph_v_2[:,2] + sph_v_3[:,2])/3.
#               
#    pent_faces, pent_ind = s.find_pentagon()
#    
#    val = 1
#    phi_avg[pent_ind] = val
#    phi_avg[phi_avg !=2] = 0
#    
#    next_tri, vconn = s.next_row_tri_v(pent_ind, pent_ind) #Defines vertical connections 
#    phi_avg[next_tri] = val + 1
#    face_cache = np.hstack((pent_ind, next_tri)) 
#    next_tri2, hconn, x_v_conn=  s.next_row_tri_h(next_tri, face_cache) 
#    vconn = vconn + x_v_conn
#    phi_avg[next_tri2] = val + 1
#    for i in range(21 + 21*(s.recursion_level-4)):
#        
#        face_cache = np.hstack((face_cache, next_tri2))
#        next_tri3, vconn2 =  s.next_row_tri_v(next_tri2, face_cache)
#        phi_avg[next_tri3] = val + i + 2
#        face_cache = np.hstack((face_cache, next_tri3))
#        next_tri4, hconn2 , x_v_conn=  s.next_row_tri_h(next_tri3, face_cache)
#        phi_avg[next_tri4] = val + i + 2
#        next_tri2 = next_tri4
#        vconn = vconn + vconn2 + x_v_conn#np.hstack((vconn, vconn2))
#        #vconn = {**vconn, **vconn2}
#        #hconn = {**hconn, **hconn2}
#        hconn = hconn + hconn2#np.hstack((hconn, hconn2))
#    
#    face_cache = np.hstack((face_cache, next_tri2))
#    next_tri3, vconn2 =  s.next_row_tri_v(next_tri2, face_cache)
#    vconn = vconn + vconn2
#    phi_avg[next_tri3] = val + 21 + 21*(s.recursion_level-4) + 2
#    colours = phi_avg#sph_v_3[:,2] #Making the colors of the faces correspond to the value of the angle
#    colours = phi_avg/sum(phi_avg)
#    #rang =range(len(phi_avg));colours = np.array(rang)/sum(rang)
#    
#    #print("vconn", vconn)
#    #print("hconn", hconn)
#    s.plot_sphere(colours)   
#    



class Create_Sphere_and_Connections:
    
    
    def __init__(self):
        
        
        
        self.t = (1.0 + (5.0)**(0.5) )/ 2.0
        
        self.icosahedron_vertices = [[-1,  self.t,  0],
                                    [ 1,  self.t,  0],
                                    [-1, -self.t,  0],
                                    [1, -self.t,  0],
                                    [ 0, -1,  self.t],
                                    [0,  1,  self.t],
                                    [ 0, -1, -self.t],
                                    [ 0,  1, -self.t],
                                    [ self.t,  0, -1],
                                    [self.t,  0,  1],
                                    [-self.t,  0, -1],
                                    [-self.t,  0,  1]]
    
        self.faces = np.array(
                    [[0, 11, 5],
                    [0, 5, 1],
                    [0, 1, 7],
                    [0, 7, 10],
                    [0, 10, 11],
                    [1, 5, 9],
                    [5, 11, 4],
                    [11, 10, 2],
                    [10, 7, 6],
                    [7, 1, 8],
                    [3, 9, 4],
                    [3, 4, 2],
                    [3, 2, 6],
                    [3, 6, 8],
                    [3, 8, 9],
                    [4, 9, 5],
                    [2, 4, 11],
                    [6, 2, 10],
                    [8, 6, 7],
                    [9, 8, 1]])
        
        self.s = sp.Sphere(vertices = self.icosahedron_vertices, faces = self.faces, recursion_level = 6 )
        self.s.ch = self.s.construct_icosphere()
        x, y, z   = self.s.ch.points[self.s.ch.vertices][:,0],self.s.ch.points[self.s.ch.vertices][:,1], self.s.ch.points[self.s.ch.vertices][:,2]
        vertex1, vertex2, vertex3 = self.s.ch.points[self.s.ch.simplices[:,0]],self.s.ch.points[self.s.ch.simplices[:,1]], self.s.ch.points[self.s.ch.simplices[:,2]]# hull_pts[ch.simplices[:,1]], hull_pts[ch.simplices[:,2]]
        
        
        
        #CHANGE TO SPHERICAL POLARS
        sph_v_1 = self.s.sph_polar_convert(vertex1)
        sph_v_2 = self.s.sph_polar_convert(vertex2) 
        sph_v_3 = self.s.sph_polar_convert(vertex3) 
        
        phi_avg = (sph_v_1[:,2] + sph_v_2[:,2] + sph_v_3[:,2])/3.
        
        #s.plot_sphere(phi_avg)
                
        pent_faces, self.pent_ind = self.s.find_pentagon()
        
        val = 1
        phi_avg[self.pent_ind] = val
        phi_avg[phi_avg !=2] = 0
        
        next_tri, vconn = self.s.next_row_tri_v(self.pent_ind, self.pent_ind) #Defines vertical connections 
        
        phi_avg[next_tri] = val + 1
        face_cache = np.hstack((self.pent_ind, next_tri)) 
        next_tri2, hconn, x_v_conn=  self.s.next_row_tri_h(next_tri, face_cache) 
        vconn = vconn + x_v_conn
        phi_avg[next_tri2] = val + 1
        
        for i in range(21*(self.s.recursion_level-3)+ 3):
            
            face_cache = np.hstack((face_cache, next_tri2))
            next_tri3, vconn2 =  self.s.next_row_tri_v(next_tri2, face_cache)
            phi_avg[next_tri3] = val + i + 2
            face_cache = np.hstack((face_cache, next_tri3))
            next_tri4, hconn2 , x_v_conn=  self.s.next_row_tri_h(next_tri3, face_cache)
            phi_avg[next_tri4] = val + i + 2
            next_tri2 = next_tri4
            vconn = vconn + vconn2 + x_v_conn#np.hstack((vconn, vconn2))
            #vconn = {**vconn, **vconn2}
            #hconn = {**hconn, **hconn2}
            self.hconn = hconn + hconn2#np.hstack((hconn, hconn2))
        
        face_cache = np.hstack((face_cache, next_tri2))
        next_tri3, vconn2 =  self.s.next_row_tri_v(next_tri2, face_cache)
        self.vconn = vconn + vconn2
        phi_avg[next_tri3] = val + 7*(self.s.recursion_level-1) + 2
        
        self.colours = phi_avg/sum(phi_avg)


#####################################################################

h= Create_Sphere_and_Connections()
n = create_network(array_nodesindices = np.arange(len(h.colours)),
                   array_vertical = h.vconn,
                   array_transv = h.hconn,
                   p_transv = 1,
                   impulse_start = h.pent_ind,
                   p_dysf = 0.8,
                   p_unexcitable = 0.05)
runc = run(network = n,
           plot=True,store=False,runs=30)
runc.animator(h.s)
#"""
#colours = n.nodes
#colours = colours/sum(colours)
#s.plot_sphere(colours)


    
    
    
    