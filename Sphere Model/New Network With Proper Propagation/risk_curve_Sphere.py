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
import Network_Projection_appended as Network

kishmatrix = np.array([[0.02,0.99981,4.3015e-06],
[0.04,0.99983,3.8088e-06],
[0.06,0.9998,1.0454e-05],
[0.08,0.99968,3.0663e-05],
[0.11,0.99772,0.00044859],
[0.13,0.96099,0.018246],
[0.15,0.60984,0.054379],
[0.17,0.16381,0.041092],
[0.19,0.017807,0.0080603],
[0.21,0.020737,0.016513],
[0.23,4.922e-05,4.8685e-05],
[0.25,0.0001084,8.4968e-05],
[0.27,0,0],
[0.29,0,0],
[0.12,0.99152,0.0027053],
[0.14,0.86184,0.028043],
[0.16,0.29714,0.055185],
[0.18,0.039206,0.013863],
[0.2,0.0056277,0.0028284],
[0.22,4.834e-05,3.6005e-05],
[0.24,0.00082172,0.00081342],
[0.26,0,0],
[0.28,0,0],
[0.3,9.406e-05,9.3115e-05],
[0.1,0.99919,0.00010423]])
kishnu = kishmatrix[:,0]
p_transversalconn= kishnu



timet=10000
number_of_systems=50
prob_dysf=0.05
prob_unexcitable=0.05
excitationvalue=25
heart=110
threshold=350



#opening pickled files for reinstatement of connections
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

conn=[]  #stored lists of connections for esch system
dysfarray=[] #stored lists of dysf array for each system
num_cells=[] #number of excited cells each time   # WATCH OUT!!!!!!!!!!!!!!!!!!!!!!!!!!!a lot of memory wasted, comment for not trial

 
t_fib_data = []#times in fibrillation data
t_in_fib = [] #actual time in fibrillation for each system


risk=[]      #risk final output
riskerror=[]   #risk std

for elements in p_transversalconn:
    n = Network.create_network(array_nodesindices = np.arange(len(colours)),
                   array_vertical = vconn,
                   array_transv = hconn,
                   p_transv = p_transversalconn,
                   impulse_start = pent_ind,
                   p_dysf = prob_dysf,
                   p_unexcitable = prob_unexcitable,
                   excitation = excitationvalue, 
                   hbs = heart)
    
    
    
    for i in range(number_of_systems):
        
        
        
        
        
        runc = Network.run(network = n, plot=False,store=True,runs=timet,fib_threshold=threshold)
        conn.append(n.connections)
        dysfarray.append(n.dysf)
        runc.propagate_storage()
        num_cells.append(runc.num_excited)  # WATCH OUT!!!!!!!!!!!!!!!!!!!!!!!!!!!a lot of memory wasted, comment for not trial
        
        
        #global storage variables to be  saved
        t_fib_data.append( runc.tfibrillation)  
        t_in_fib.append(runc.timeinfibrillation())
        
        
        n.reinitialise()
    
    risk.append(np.average(  t_in_fib[-number_of_systems:] ) ) #appends performs the average on the last number of systems elements in the t_fib_data list which correspond to one value of the risk curve
    riskerror.append(np.std(  t_in_fib[-number_of_systems:] )/ (number_of_systems**0.5))



num=open('excited_cells', 'wb')   # WATCH OUT NEXT THREE LINES!!!!!!!!!!!!!!!!!!!!!!!!!!!a lot of memory wasted, comment for not trial
pickle.dump(num_cells, num, -1) 
num.close()



fileobj0 = open('connections.pkl', 'wb')
fileobj1 = open('dysf_grids.pkl', 'wb')
fileobj2 = open('t_fib_data.pkl', 'wb')
fileobj3 = open('t_in_fib.pkl', 'wb')
fileobj4 = open('risk.pkl', 'wb')
fileobj5 = open('riskerror.pkl', 'wb')


pickle.dump(conn, fileobj0, -1)
pickle.dump(dysfarray, fileobj1, -1)
pickle.dump(t_fib_data, fileobj2, -1)
pickle.dump(t_in_fib, fileobj3, -1)
pickle.dump(risk, fileobj4, -1)
pickle.dump(riskerror, fileobj5, -1)

 
fileobj0.close()
fileobj1.close()
fileobj2.close()
fileobj3.close()
fileobj4.close()
fileobj5.close()


plt.figure()
plt.xlabel("Percentage transversal Connections/Nu")
plt.ylabel("Mean time in AF/Risk of AF")
plt.errorbar(p_transversalconn,risk, yerr=riskerror, fmt='o')
plt.title('Risk Curve')
plt.show()           
        

