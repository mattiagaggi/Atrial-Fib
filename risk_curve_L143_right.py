import FinalAF_L143_right as FinalAF
import numpy as np
from matplotlib import pyplot as plt
import pickle

#Kishan'f Data
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
#some_nu = [0.31, 0.32, 0.34,0.33, 0.35, 0.36, 0.37, 0.38,0.39, 0.4, 0.41, 0.42, 0.43,0.44, 0.45]

extra_nu = [0.195, 0.205,0.215, 0.225, 0.235, 0.245,0.32, 0.34,0.33, 0.35, 0.36, 0.37, 0.38,0.39, 0.4]
m_trans_conn = p_transversalconn.copy()
for nu in extra_nu:
    m_trans_conn.append(nu)

#m_trans_conn = extra_nu
time=100000
number_of_systems=50




dysfgrids=[]
fibrosisgrids=[]

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

for elements in m_trans_conn:
    print("Nu = %s" %(elements))

    average=[]
    singlemeasures=[]
    t_fib_temp = [] #list of times in fibrillation
    af_number = []  #list of fibrillation times
    number_fib = [] #number of systems in fibrillation
  
    systems=FinalAF.heart(143,0.05,elements,0.05)

   
    for i in range(number_of_systems):
        print(i)
        
        dysfgrids.append(systems.dysfgrid)
        fibrosisgrids.append(systems.edgegrid)
    
        systemsrun=FinalAF.run(systems,False,True,time,False)
        
        timeinfib = systemsrun.timeinfibrillation()
        listoffib = systemsrun.tfibrillation
        n_fib = len(listoffib)
        systems.reinitialise()
        t_fib_temp.append(timeinfib)
        

        singlemeasures.append(timeinfib/float(time))
        af_number.append(listoffib)
        if n_fib > 0:
            number_fib.append(1)
        else:
            number_fib.append(0)
     
    t_in_fib.append(t_fib_temp)
    t_fib_data.append(af_number)
    
    n_fib_data.append(number_fib)
    n_fib_avg.append(np.mean(np.array(number_fib)))
    n_fib_err.append(np.std(np.array(number_fib))/((number_of_systems)**(0.5)))
    average=np.mean(np.array(singlemeasures))
    std=np.std(np.array(singlemeasures))/((len(singlemeasures))**(0.5)) #Standard Error of the mean
    risk.append(average)
    riskstd.append(std)

risk=np.array(risk)

fileobj0 = open('listoffibrillationtimes143_right.pkl', 'wb')
fileobj1 = open('risk143_right.pkl', 'wb')
fileobj2 = open('riskerrorr143_right.pkl', 'wb')
fileobj4 = open('numberfibrillations143_right.pkl', 'wb')
fileobj5 = open('nfibaverage143_right.pkl', 'wb')
fileobj6 = open('nfiberr143_right.pkl', 'wb')
fileobj7 = open('tinfib143_right.pkl', 'wb')
fileobj8 = open('dysfgrids143_right.pkl', 'wb')
fileobj9 = open('fibrosisgrid143s_right.pkl', 'wb')

pickle.dump(t_fib_data, fileobj0, -1)
pickle.dump(risk, fileobj1, -1)
pickle.dump(riskstd, fileobj2, -1)
pickle.dump(n_fib_data, fileobj4, -1)
pickle.dump(n_fib_avg, fileobj5, -1)
pickle.dump(n_fib_err, fileobj6, -1)
pickle.dump(t_in_fib, fileobj7, -1)
pickle.dump(dysfgrids, fileobj8, -1)
pickle.dump(fibrosisgrids, fileobj9, -1)

fileobj1.close()
fileobj2.close()
fileobj4.close()
fileobj0.close()
fileobj5.close()
fileobj6.close()
fileobj7.close()
fileobj8.close()
fileobj9.close()
    
#plt.figure()
def p_risk_sph_f(nu, tau, rec, delta):
    val = 1- (1-(1-nu)**(tau))**( 5*delta*(4**(rec)))# -640 )  - (1-(1-nu)**((tau/2)*(3)))**(640*delta) )
    return val

def p_risk(nu, L, tau, delta):
    val = 1- ( 1-(1-nu)**(tau) )**( delta*L*L)# -640 )  - (1-(1-nu)**((tau/2)*(3)))**(640*delta) )
    return val

import numpy as np
nulist = np.linspace(0,0.45, 100)
tau = 25
tau2 = 35
rec = 6   
delta = 0.05
L=143
an_p_risk = []
p_risk_n = []
for i in nulist:
    #i = 1-(1-i)**2
    prisk1 = p_risk_sph_f(i, tau,rec, delta)
    prisk2 = p_risk(i, L, tau2, delta)
    #print(prisk)
    an_p_risk.append(prisk1)
    p_risk_n.append(prisk2)
    
fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel("Percentage transversal Connections/Nu")
ax.set_ylabel("Mean time in AF/Risk of AF")

ax.plot(m_trans_conn,risk, 'b+', label = 'Original Model - L=143')
ax.errorbar(m_trans_conn,risk, yerr=riskstd, fmt='b+')

ax.plot(kishnu, kishrisk, 'g^', label = 'Manani/Christensen Data')
ax.errorbar(kishnu, kishrisk, yerr=kisherror, fmt='g^')

#plt.xlabel("Percentage transversal Connections/Nu")
#plt.ylabel("Mean time in AF/Risk of AF")
#plt.plot(p_transversalconn,risk, 'bo')
#plt.errorbar(m_trans_conn,risk, yerr=riskstd, fmt='o')
#plt.errorbar(p_transversalconn,n_fib_avg, yerr=n_fib_err, fmt='o')
#plt.plot(kishnu, kishrisk, 'ro')
ax.plot(nulist, p_risk_n, 'r--', label='Analytic Risk Curve Christensen')
ax.plot(nulist, an_p_risk, 'r--', label='Analytic Risk Curve Sphere')


ax.legend()
ax.set_title('Risk Curve--Christensen Model')
plt.show()    
 
        
        
        