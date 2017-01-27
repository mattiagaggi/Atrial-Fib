import numpy as np
from scipy.sparse import csr_matrix
import operator




def decision(p):
    if np.random.random()<p:
        return True
    else:
        return False

class create_network:
    def __init__(self,array_nodesindices,array_vertical,array_transv,p_transv,p_dysf,p_unexcitable):
        """1- create a dictionary from the tuples with repeated elements with indices - DONE
            2-define excited cells in an array or list of indices - DONE
            3- check that the dictionary does the proper mapping --- what is the quickest way to do the mapping?
            4 -check if the cell is dysfunctional
            5 check if the cell is refractory"""
            
        
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
     
    def excite(self,nodeindex):
        
        self.node[nodeindex]+=self.excitation
        self.excited.append(nodeindex)

    def onestep(self):
        
        newnodes=np.zeros(self.size)
        f = operator.itemgetter(*self.excited) # works well with list or array
        self.excited=f(self.connections) #the output of f() is ( [....],[...],....) change into a unique list
        newnodes[self.excited]+=1
        # alternative function: list(self.connections[_] for _ in self.excited) check which one is faster
        
        self.unexcited=np.random.rand(self.size)  #create array of dysf cells maybe can find a quicker way
        self.unexcited[self.unexcited < self.p_unexcitable] = 1
        self.unexcited[self.self.unexcited != 1] = 0
        
        newnodes-=self.dysf*self.unexcited*newnodes #remove excited dysf cells which  are not excited
        
        newnodes*= (self.nodes==0) #removes refractory cells
        
        self.nodes-=1  
        self.nodes[self.nodes==-1]=0
        
        self.nodes += newnodes*self.excitation
   


class propagation:
    
    def __init__(self,network, runs=10):
        
        
        self.runs=runs    
        self.network=network

       
        for i in range(runs):
            
            new_nodes=np.array(self.network.adjacency*self.network.nodes)
            self.excited=np.random.rand(self.network.size,1)
            #need to erase the dysf cells from here
            self.excited*=self.dysf
            #self.excited[.....  set to 1
            new_nodes #readd cells that are dysf and excited
            
            new_nodes*=(self.network.nodes==0)#check if nodes are refractory