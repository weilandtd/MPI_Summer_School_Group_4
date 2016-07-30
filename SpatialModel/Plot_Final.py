import numpy as np 
import matplotlib.pyplot as plt

 
dt = 0.1
t = 10.0
Nt = int(t/dt)+1

folder = "PAR"
filename = "time"
Timeseries = []

color = ['b', 'g', 'r','c', 'm', 'k']
cc = 0
for i in range(0,Nt,5):
    try:
        Data = np.loadtxt(folder+"/"+filename+"_%d.txt"%(i))
        sim1, = plt.plot(Data[:,0], Data[:,2],color='blue', linewidth = 1.0, label ='Simulation' )

        sim2, = plt.plot(Data[:,0], Data[:,3],color='red', linewidth = 1.0, label ='Simulation' )
    # Plot your Data 
        #print 'Sim'
        #print np.var(Data[:,2])
    #Timeseries.append(Data)
    except: pass
    cc += 1
#Timeseries = np.array(Timeseries)
#ana = mpatches.Patch(color='red', label='')
#sim = mpatches.Patch(color='red', label='Simulation')

plt.legend([sim1, sim2], [ 'aPAR','pPAR' ])
plt.ylabel('Concentration')
plt.xlabel('x-coordinate')
plt.show()


