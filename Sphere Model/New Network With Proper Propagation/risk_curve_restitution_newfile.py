import Network_Projection_withrestitution as Network
import numpy as np
from matplotlib import pyplot as plt
import pickle
import time
#Kishan's Datae
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
kishrisk = kishmatrix[:,1]
kisherror = kishmatrix[:,2]
p_transversalconn= kishnu #np.arange(0,0.6,0.05) #nu
p_transversalconn= list(kishnu) #np.arange(0,0.6,0.05) #nu
p_transversalconn.append(0.01)
p_transversalconn.append(0.005)
np.asarray(p_transversalconn)

extra_nu = [0.195, 0.205,0.215, 0.225, 0.235, 0.245,0.32, 0.34,0.33, 0.35, 0.36, 0.37, 0.38,0.39, 0.4]
m_trans_conn = p_transversalconn.copy()
for nu in extra_nu:
    m_trans_conn.append(nu)

timet=100000
number_of_systems=12
prob_dysf=0.05
prob_unexcitable=0.05
excitationvalue=25
heart=110
threshold=350



#fileobj0 = open('realizationdata', 'wb')

d=open('dysf_grids_sph_run1.pkl','rb')
dysfgrids=pickle.load(d)
d.close()

dx=open('dysf_grids_sph_run1_extra_data.pkl','rb')
dysfgrids_x=pickle.load(dx)
dx.close()

c=open('connections_sph_run1.pkl','rb')
connections=pickle.load(c)
c.close()

cx=open('connections_sph_run1_extra_data.pkl','rb')
connections_x=pickle.load(cx)
cx.close()

#f=open('fibrosisgrids.pkl','rb')
#fibrosisgrids=pickle.load(f)
#f.close()



counter=0

"""
if len(dysfgrids)!=len(p_transversalconn)*len(number_of_systems):
    "watch it mate dysfgrids has wrong number  of elements maybe"
if len(fibrosisgrids)!=len(p_transversalconn)*len(number_of_systems):
    "watch it mate dysfgrids has wrong number  of elements maybe"
    
"""

d=open('horiz_conn_rec_6xconn.pkl', 'rb')
hconn=pickle.load(d)
f=open('vert_conn_rec_6xconn.pkl', 'rb')
vconn=pickle.load(f)
e=open('startimp_ind_rec_6xconn.pkl', 'rb')
pent_ind=pickle.load(e)
g=open('colours_rec_6xconn.pkl', 'rb')
colours=pickle.load(g)

d.close()
e.close()
f.close()
g.close()

for dysfg in dysfgrids_x:
    dysfgrids.append(dysfg)

for i in connections_x:
    connections.append(i)
#h=open('sph_rec_6xconn.pkl', 'rb')
#s=pickle.load(h)

#conn=[]  #stored lists of connections for esch system
#dysfarray=[] #stored lists of dysf array for each system
#num_cells=[] #number of excited cells each time   # WATCH OUT!!!!!!!!!!!!!!!!!!!!!!!!!!!a lot of memory wasted, comment for not trial

 #if run for the first time
 
 
 
 
 
 
 #cretes a variable that keeps count of the number of systems already analysed
systemsstudied=0
o = open('systemsstudied.pkl', 'wb')
pickle.dump(systemsstudied, o, -1)
o.close()


#for every fiuture runs - updates the variable


o = open('systemsstudied.pkl', 'rb')
start=pickle.load(o)
o.close()
end=start+number_of_systems
o = open('systemsstudied.pkl', 'wb')
pickle.dump(end, o, -1)
o.close()





for j in range(len(m_trans_conn)):
    
    #only when run for the first time it creates the files
    t_fib_data= []
    t_in_fib=[]
    name1='t_fib_data_rest_transversalconnnumber '+str(j)+'.pkl'
    name2='t_in_fib_sph_rest_transversalconnnumber '+str(j)+'.pkl'
    fileobj1 = open(name1, 'wb')
    fileobj2 = open(name2, 'wb')
    pickle.dump(t_fib_data, fileobj1, -1)
    pickle.dump(t_in_fib, fileobj2, -1)
    fileobj1.close()
    fileobj2.close()
    
    
    
    #for every future runs, it opens the file,reads them
    name1='t_fib_data_rest_transversalconnnumber '+str(j)+'.pkl'
    name2='t_in_fib_sph_rest_transversalconnnumber '+str(j)+'.pkl'
    fileobj1 = open(name1, 'rb')
    fileobj2 = open(name2, 'rb')
    t_fib_data= pickle.load(fileobj1)
    t_in_fib=pickle.load(fileobj2)
    fileobj1.close()
    fileobj2.close()
    
   
    
   

    n = Network.create_network(array_nodesindices = np.arange(len(colours)),
                   array_vertical = vconn,
                   array_transv = hconn,
                   p_transv = m_trans_conn[j],
                   impulse_start = pent_ind,
                   p_dysf = prob_dysf,
                   p_unexcitable = prob_unexcitable,
                   excitation = excitationvalue, 
                   hbs = heart,
                   recursion = 6)

    
    
    print("nu = ",m_trans_conn[j])
    
    for i in range(start,end):
        
        n.connections = connections[i + j*50]
        n.dysf = dysfgrids[i + j*50]
        print("system number", i)
        timea = time.time()
        
        runc = Network.run(network = n, plot=False,store=True,runs=timet,fib_threshold=threshold)
        
        #dysfarray.append(n.dysf)
        runc.propagate_storage()
        #num_cells.append(runc.num_excited)  # WATCH OUT!!!!!!!!!!!!!!!!!!!!!!!!!!!a lot of memory wasted, comment for not trial
        
        
        #global storage variables to be  saved
        t_fib_data.append( runc.tfibrillation)  
        tinfib = runc.timeinfibrillation()
        t_in_fib.append(tinfib)
        
     
        
        n.reinitialise()
        timeb = time.time()
        print("runtime for this system = " ,(timeb-timea))
    
    
    
    #it overwrite the files
    fileobj1 = open(name1, 'wb')
    fileobj2 = open(name2, 'wb')
    pickle.dump(t_fib_data, fileobj1, -1)
    pickle.dump(t_in_fib, fileobj2, -1)
    fileobj1.close()
    fileobj2.close()
       
   
    
    
    
#risk.append(np.average( t_in_fib[-number_of_systems:]/float(timet) ) ) #appends performs the average on the last number of systems elements in the t_fib_data list which correspond to one value of the risk curve
#riskerror.append(np.std( t_in_fib[-number_of_systems:]/float(timet) )/ (number_of_systems**0.5))


#nfibavg = 
#risk=np.array(risk)



        
        
        
        