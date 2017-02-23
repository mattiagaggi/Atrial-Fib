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

risk=[]
riskstd=[]
t_fib_data = []
t_fib_temp = []
number_fib = []
n_fib_data = []
n_fib_avg = []
n_fib_err = []
t_in_fib = []

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

er=open('risk.pkl','rb')
risk_o=pickle.load(er)
er.close()

fr=open('riskerrorr.pkl','rb')
risk_o_std=pickle.load(fr)
fr.close()
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

 
t_fib_data = []#times in fibrillation data
t_in_fib = [] #actual time in fibrillation for each system
risk_fib = []

risk=[]      #risk final output
riskerror=[]   #risk std


for j in range(len(m_trans_conn)):
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
    for i in range(number_of_systems):
        
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
        risk_fib.append(tinfib/float(timet))
        
        
        n.reinitialise()
        timeb = time.time()
        print("runtime for this system = " ,(timeb-timea))
    
    risk.append(np.average( risk_fib[-number_of_systems:] ) ) #appends performs the average on the last number of systems elements in the t_fib_data list which correspond to one value of the risk curve
    riskerror.append(np.std( risk_fib[-number_of_systems:] )/ (number_of_systems**0.5))


#nfibavg = 
risk=np.array(risk)



fileobj2 = open('t_fib_data_sph_rest_run_12iter.pkl', 'wb')
fileobj3 = open('t_in_fib_sph_rest_run_12iter.pkl', 'wb')
fileobj6 = open('risk_fib_sphrest_run_12iter.pkl', 'wb')
fileobj4 = open('risk_sph_rest_run_12iter.pkl', 'wb')
fileobj5 = open('riskerror_sph_rest_run_12iter.pkl', 'wb')



pickle.dump(t_fib_data, fileobj2, -1)
pickle.dump(t_in_fib, fileobj3, -1)
pickle.dump(risk, fileobj4, -1)
pickle.dump(riskerror, fileobj5, -1)
pickle.dump(risk_fib, fileobj6, -1)
 

fileobj2.close()
fileobj3.close()
fileobj4.close()
fileobj5.close()
fileobj6.close()




file1 = open('risk.pkl', 'rb')
risk_o = pickle.load(file1)
file1.close()

file2 = open('riskerrorr.pkl', 'rb')
risk_o_std = pickle.load(file2)
file2.close()

file3 = open('risk_restn1.pkl', 'rb')
risk_rest = pickle.load(file3)
file3.close()

file4 = open('riskerrorr_restn1.pkl', 'rb')
riskstd_rest = pickle.load(file4)
file4.close()


file5 = open('risk_sph_run1.pkl', 'rb')
risk_sph_run1 = pickle.load(file5)
file5.close()

file6 = open('risk_sph_run1_extra_data.pkl', 'rb')
riskx_sph_run1 = pickle.load(file6)
file6.close()

file7 = open('riskerror_sph_run1.pkl', 'rb')
riskstd_sph_run1 = pickle.load(file7)
file7.close()

file8 = open('riskerror_sph_run1_extra_data.pkl', 'rb')
riskstdx_sph_run1 = pickle.load(file8)
file8.close()

for i in riskx_sph_run1:
    risk_sph_run1.append(i)
for j in riskstdx_sph_run1:
    riskstd_sph_run1.append(j)


"""
plt.figure()
plt.xlabel("Percentage transversal Connections/Nu")
plt.ylabel("Mean time in AF/Risk of AF")
plt.errorbar(p_transversalconn,risk, yerr=riskerror, fmt='o')
plt.title('Risk Curve')
plt.show()           
"""
      
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("Percentage transversal Connections/Nu")
ax.set_ylabel("Mean time in AF/Risk of AF")

ax.plot(p_transversalconn,risk_rest, 'bo', label = 'Restitution Model')
ax.errorbar(p_transversalconn,risk_rest, yerr=riskstd_rest, fmt='bo')

ax.plot(p_transversalconn,risk_o, 'r+', label = 'Original Model')
ax.errorbar(p_transversalconn,risk_o, yerr=risk_o_std, fmt='r+')
#plt.errorbar(p_transversalconn,n_fib_avg, yerr=n_fib_err, fmt='o')
ax.plot(kishnu, kishrisk, 'g^', label = 'Kishans Data')
ax.errorbar(kishnu, kishrisk, yerr=kisherror, fmt='g^')

ax.plot(m_trans_conn, risk_sph_run1, 'cd', label = 'Sphere Model')
ax.errorbar(m_trans_conn,risk_sph_run1, yerr=riskstd_sph_run1, fmt='cd')

ax.plot(m_trans_conn,risk, 'm*', label = 'Sphere Model with Restitution')
ax.errorbar(m_trans_conn,risk, yerr=riskerror, fmt='m*')


ax.legend()
ax.set_title('Risk Curve--Sphere Model with Restitution')
plt.show()    
        
        
        
        