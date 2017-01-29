import numpy as np
from scipy.sparse import csr_matrix
import operator
from itertools import chain
import time


def decision(p):
    if np.random.random()<p:
        return True
    else:
        return False


class create_network:
    def __init__(self,array_nodesindices,array_vertical,array_transv,p_transv,p_dysf,p_unexcitable):
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
    
        self.array_vertical=array_vertical
        self.size=len(array_nodesindices)
        self.nodes=np.zeros(self.size)
       
        self.array_transv=[]
        self.excited=[]

        self.dysf=np.random.rand(self.size)  #create array of dysf cells
        self.dysf[self.dysf < self.p_dysf] = 1
        self.dysf[self.dysf != 1] = 0
        self.unexcited=np.random.rand(self.size)
        self.array_alltransv=array_transv
               
        for elements in array_transv: #append transversal connections to a new list according to the probability of transv connect
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
        
        self.nodes[nodeindex]+=self.excitation
        self.excited.append(nodeindex)

    def onestep(self):
        """
        oper 1-2: maps excited elem to the next excited
        
        oper 3: excites the newnodes
        
        oper 4 : creates array of unexcited elements
        
        oper5: removes nodes in refractory state, then removes dysfunctional non-excited elements
        
        oper 6: every node-=1 and adds the newnodes
            
        """ 
        
        
        time0=time.time()
        newnodes=np.zeros(self.size)
        time1=time.time()
    
        
        
        
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
        newnodes[newnodes>0]=1#need to account for repeated elements
        time4=time.time()
        print("time taken",(time4-time3))
    
        self.unexcited=np.random.rand(self.size)  #create array of dysf cells maybe can find a quicker way #this seems to take long 0.001
        self.unexcited[self.unexcited < self.p_unexcitable] = 1
        self.unexcited[self.unexcited != 1] = 0
        time5=time.time()
        print("time taken",(time5-time4))
        
        newnodes=newnodes-(self.dysf*self.unexcited*newnodes) #remove excited dysf cells which  are not excited
        newnodes*= (self.nodes==0) #removes refractory cells
        time6=time.time()
        print("time taken",(time6-time5))
        
        self.nodes-=1  
        self.nodes[self.nodes==-1]=0
        self.nodes += newnodes*self.excitation
        time7=time.time()
        print("time taken",(time7-time6))
        
        print ("total time", ( (time2-time1)+(time3-time2)+(time4-time3)+(time5-time4)+(time6-time5)+(time7-time6)))
   


class propagation:
    
    def __init__(self,network, runs=10):
        
        
        self.runs=runs    
        self.network=network

       
        
            
            
        