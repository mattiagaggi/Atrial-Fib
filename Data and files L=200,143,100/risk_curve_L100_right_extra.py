import FinalAF_L100_right as FinalAF
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

some_new_nu = [0.41,0.42, 0.43, 0.44, 0.45]
zeross = [0,0,0,0,0]
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

for elements in some_new_nu:
    print("Nu = %s" %(elements))

    average=[]
    singlemeasures=[]
    t_fib_temp = [] #list of times in fibrillation
    af_number = []  #list of fibrillation times
    number_fib = [] #number of systems in fibrillation
  
    systems=FinalAF.heart(100,0.05,elements,0.05)

   
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

fileobj0 = open('listoffibrillationtimes_L100_extra.pkl', 'wb')
fileobj1 = open('risk_L100_extra.pkl', 'wb')
fileobj2 = open('riskerrorr_L100_extra.pkl', 'wb')
fileobj4 = open('numberfibrillations_L100_extra.pkl', 'wb')
fileobj5 = open('nfibaverage_L100_extra.pkl', 'wb')
fileobj6 = open('nfiberr_L100_extra.pkl', 'wb')
fileobj7 = open('tinfib_L100_extra.pkl', 'wb')
fileobj8 = open('dysfgrids_L100_extra.pkl', 'wb')
fileobj9 = open('fibrosisgrid_L100_extra.pkl', 'wb')

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

nulist = np.linspace(0,0.45, 100)
tau = 25
tau2 = 35
tau3 = 50
rec = 6   
delta = 0.05
L=143
L2 = 100
L3 = 200
an_p_risk = []
p_risk_n_143 = []
p_risk_n_100 = []
p_risk_n_200 = []
for i in nulist:
    #i = 1-(1-i)**2
    prisk1 = p_risk_sph_f(i, tau,rec, delta)
    prisk2 = p_risk(i, L, tau2, delta)
    prisk3 = p_risk(i, L2, tau, delta)
    prisk4 = p_risk(i, L3, tau3, delta)
    #print(prisk)
    an_p_risk.append(prisk1)
    p_risk_n_143.append(prisk2)
    p_risk_n_100.append(prisk3)
    p_risk_n_200.append(prisk4)
    

#mt2 = mt[mt<=0.22]
#for j in range(len(mt2)):

#for k in m_trans_conn2:
    
#risk
m_trans_conn2 = np.asarray(m_trans_conn)
m_trans_conn2 = m_trans_conn2[m_trans_conn2 >0.22]

  
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



fileobj1 = open('risk_restn1_143_trial2.pkl', 'rb')
risk_restn1_143 = pickle.load(fileobj1)
fileobj1.close()

fileobj2 = open('riskerrorr_restn1_143_trial2.pkl', 'rb')
risk_restn1_143_error = pickle.load(fileobj2)
fileobj2.close()

fileobj1 = open('risk143_right.pkl', 'rb')
risk_o_143 = pickle.load(fileobj1)
fileobj1.close()

fileobj2 = open('riskerrorr143_right.pkl', 'rb')
risk_o_std_143 = pickle.load(fileobj2)
fileobj2.close()

fileobj1 = open('risk_L100.pkl', 'rb')
risk_L100 = pickle.load(fileobj1)
fileobj1.close()

fileobj2 = open('riskerrorr_L100.pkl', 'rb')
risk_o_std_100 = pickle.load(fileobj2)
fileobj2.close()



fileobj2 = open('nu_sphere_restitution.pkl', 'rb')
nu_sphere_restitution= pickle.load(fileobj2)
fileobj2.close()


fileobj2 = open('risk_error_sphere_restitution_final.pkl', 'rb')
risk_error_sphere_restitution_final = pickle.load(fileobj2)
fileobj2.close()


fileobj2 = open('risk_sphere_restitution_final.pkl', 'rb')
risk_sphere_restitution_final= pickle.load(fileobj2)
fileobj2.close()




for i in range(len(riskstdx_sph_run1)):
    risk_sph_run1.append(riskx_sph_run1[i])
    riskstd_sph_run1.append(riskstdx_sph_run1[i])


#####################################################################################
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.set_xlabel("Percentage transversal Connections/Nu", fontsize = 20)
ax1.set_ylabel("Mean time in AF/Risk of AF", fontsize = 20)

ax1.plot(p_transversalconn,risk_o, 'rs', label = 'Original Model L=200')
ax1.errorbar(p_transversalconn,risk_o, yerr=risk_o_std, fmt='rs')
#plt.errorbar(p_transversalconn,n_fib_avg, yerr=n_fib_err, fmt='o')
ax1.plot(p_transversalconn,risk_rest, 'bo', label = 'Restitution Model L=200')
ax1.errorbar(p_transversalconn,risk_rest, yerr=riskstd_rest, fmt='bo')

ax1.plot(kishnu, kishrisk, 'g^', label = 'Manani/Christensen Data L=200')
ax1.errorbar(kishnu, kishrisk, yerr=kisherror, fmt='g^')

ax1.plot(nulist, p_risk_n_200, 'r--', label='Analytic Risk Curve L=200 Christensen')

ax1.plot(some_new_nu, zeross, 'rs')
ax1.plot(some_new_nu, zeross, 'bo')
ax1.plot(some_new_nu, zeross, 'mh')
ax1.plot(some_new_nu, zeross, 'y>')
ax1.plot(some_new_nu, zeross, 'cd')
ax1.plot(some_new_nu, zeross, 'g*')

ax1.plot(m_trans_conn,risk_o_143, 'mh', label = 'Original Model L=143')
ax1.errorbar(m_trans_conn,risk_o_143, yerr=risk_o_std_143, fmt='mh')

ax1.plot(m_trans_conn,risk_restn1_143, 'y>', label = 'Restitution Model L=143')
ax1.errorbar(m_trans_conn,risk_restn1_143, yerr=risk_restn1_143_error, fmt='y>')

ax1.plot(nulist, p_risk_n_143, 'm--', label='Analytic Risk Curve L=143 Christensen')

fileobj1 = open('risk_restn1_1_L100.pkl', 'rb')
risk_L100_rest=pickle.load(fileobj1)
fileobj1.close()

fileobj2 = open('riskerrorr_restn1__L100.pkl', 'rb')
risk_o_std_100_rest=pickle.load(fileobj2)
fileobj2.close()

ax1.plot(m_trans_conn,risk_L100, 'kh', label = 'Original Model L=100')
ax1.errorbar(m_trans_conn,risk_L100, yerr=risk_o_std_100, fmt='kh')

ax1.plot(m_trans_conn,risk_L100_rest, 'r<', label = 'Restitution Model L=100')
ax1.errorbar(m_trans_conn,risk_L100_rest, yerr=risk_o_std_100_rest, fmt='r<')

ax1.plot(some_new_nu,risk, 'kh' ) #EXTRADATA
ax1.errorbar(some_new_nu,risk, yerr=riskstd) #EXTRADATA

#ax1.plot(m_trans_conn,risk, 'r<', label = 'Restitution Model L=100') #EXTRADATA
#ax1.errorbar(m_trans_conn,risk, yerr=riskstd, fmt='r<') #EXTRADATA

ax1.plot(nulist, p_risk_n_100, 'k--', label='Analytic Risk Curve L=143 Christensen')



ax1.plot(m_trans_conn, risk_sph_run1, 'cd', label = 'Sphere Model')
ax1.errorbar(m_trans_conn,risk_sph_run1, yerr=riskstd_sph_run1, fmt='cd')

ax1.plot(nu_sphere_restitution,risk_sphere_restitution_final, 'g*', label = 'Sphere Model with Restitution')
ax1.errorbar(nu_sphere_restitution,risk_sphere_restitution_final, yerr=risk_error_sphere_restitution_final, fmt='g*')



#ax.plot(m_trans_conn,riskv2, 'm*', label = 'Sphere Model with Restitution')
#ax.errorbar(m_trans_conn,riskv2, yerr=riskerrv2, fmt='m*')

ax1.plot(nulist, an_p_risk, 'c--', label='Analytic Risk Curve Sphere')  

#ax1.legend()
#ax1.set_title('Risk Curve--All Models')
#plt.show()    

    
#fig = plt.figure()
#ax = ax1
#ax.set_xlabel("Percentage transversal Connections/Nu")
#ax.set_ylabel("Mean time in AF/Risk of AF")

#ax.plot(nulist, an_p_risk, 'r--', label='Analytic Risk Curve Sphere')
from matplotlib.font_manager import FontProperties

fontP = FontProperties()
fontP.set_size('medium')
ax1.legend()
legend = ax1.legend(loc=0, ncol=1, bbox_to_anchor=(0, 0, 1, 1),
           prop = fontP,fancybox=True,shadow=False,title='LEGEND')

plt.setp(legend.get_title(),fontsize='medium')
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
#ax.set_title('Risk Curve L=143')
plt.show()
        
        