import AF_restitution
import numpy as np
from matplotlib import pyplot as plt
import pickle

#Kishan's Data
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
time=10000
number_of_systems=25
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

for elements in p_transversalconn:
    print("Nu = %s" %(elements))
    totfibtime = 0
    average=[]
    singlemeasures=[]
    t_fib_temp = []
    af_number = []
    number_fib = []
    #n_fib_real = 0
    systems=AF_restitution.heart(200,0.05,1-elements,0.05)
    #systems.electrocardiosetup([100,100])
    for i in range(number_of_systems):
        print(i)
        
        systemsrun=AF_restitution.run(systems,False,True,time,False)
        timeinfib = systemsrun.timeinfibrillation()
        listoffib = systemsrun.tfibrillation
        n_fib = len(listoffib)
        systems.reinitialise()
        t_fib_temp.append(timeinfib)
        #pickle.dump(systemsrun.gridintime, fileobj0, -1)

        singlemeasures.append(timeinfib/float(time))
        af_number.append(listoffib)
        if n_fib > 0:
            number_fib.append(1)
        else:
            number_fib.append(0)
    print "one set of systems done"
    t_in_fib.append(t_fib_temp)
    t_fib_data.append(af_number)
    n_fib_data.append(number_fib)
    n_fib_avg.append(np.mean(np.array(number_fib)))
    n_fib_err.append(np.std(np.array(number_fib))/((number_of_systems)**(0.5)))
    average=np.mean(np.array(singlemeasures))
    std=np.std(np.array(singlemeasures))/((len(singlemeasures))**(0.5)) #Standard Error of the mean
    risk.append(average)
    riskstd.append(std)





#nfibavg = 
risk=np.array(risk)
from matplotlib import pyplot as plt
plt.figure()
plt.xlabel("Percentage transversal Connections/Nu")
plt.ylabel("Mean time in AF/Risk of AF")
#plt.plot(p_transversalconn,risk, 'bo')
plt.errorbar(p_transversalconn,risk, yerr=riskstd, fmt='o')
#plt.errorbar(p_transversalconn,n_fib_avg, yerr=n_fib_err, fmt='o')
#plt.plot(kishnu, kishrisk, 'ro')
plt.errorbar(kishnu, kishrisk, yerr=kisherror, fmt='o')
plt.legend('Data from Model', 'Kishans Data')
plt.title('Risk Curve')
plt.show()   

fileobj7 = open('tinfibdatarightt2.pkl', 'wb')
fileobj0 = open('tfibdatarightt2.pkl', 'wb')
fileobj1 = open('riskrcdatarightt2.pkl', 'wb')
fileobj2 = open('riskerrorrcdatarightt2.pkl', 'wb')
fileobj3 = open('tfibrightt2.pkl', 'wb')
fileobj4 = open('nfibrightt2.pkl', 'wb')
fileobj5 = open('nfibavgrightt2.pkl', 'wb')
fileobj6 = open('nfiberrrightt2.pkl', 'wb')

pickle.dump(t_fib_data, fileobj0, -1)
pickle.dump(risk, fileobj1, -1)
pickle.dump(riskstd, fileobj2, -1)
pickle.dump(t_fib_data, fileobj3, -1)
pickle.dump(n_fib_data, fileobj4, -1)
pickle.dump(n_fib_avg, fileobj5, -1)
pickle.dump(n_fib_err, fileobj6, -1)
pickle.dump(t_in_fib, fileobj7, -1)
fileobj1.close()
fileobj2.close()
fileobj3.close()
fileobj4.close()
fileobj0.close()
fileobj5.close()
fileobj6.close()
    
        
        
        
        