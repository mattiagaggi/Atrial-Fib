import FinalAF
import numpy as np
from matplotlib import pyplot as plt
p_transversalconn=np.arange(0,0.6,0.05)
time=10000
number_of_systems=2
risk=[]
riskstd=[]

for elements in p_transversalconn:
    average=[]
    singlemeasures=[]
    for systems in range(number_of_systems):
        systems=FinalAF.heart(200,0.05,1-elements,0.05)
        systemsrun=FinalAF.run(systems,False,True,time,False)
        singlemeasures.append(systemsrun.timeinfibrillation())
    average=np.mean(np.array(singlemeasures))
    std=np.std(np.array(singlemeasures))
    risk.append(average)
    riskstd.append(std)


risk=np.array(risk)/float(time)
plt.figure()
plt.xlabel("Percentage transversal Connections")
plt.ylabel("Risk")
pl,=plt.plot(p_transversalconn,risk)
eb = plt.errorbar(p_transversalconn,risk, yerr=riskstd, color='b')
plt.show()       
        
        
        
        