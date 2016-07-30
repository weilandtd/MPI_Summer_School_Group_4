import numpy as np 
import matplotlib.pyplot as plt

#Inital distribution 
def IC(x,t):
    t0 = 0.001
    D = 1.0
    return 1.0/np.sqrt(4*np.pi*D*t0)*np.exp(-x**2/(4*D*t0))+0.1

#Solution for the diffusion
def g(x,t):
    t0 = 0.001
    D = 1.0
    return 1.0/np.sqrt(4*np.pi*D*(t0+t))*np.exp(-x**2/(4*D*(t0+t)))


def G(x,t):
    N  = 100
    G_ = 0.0
    L  = 2.0
    for n in range(-N,N+1,1):
        G_ += g(x-n*L,t)
    return G_ + 0.1

dt = 0.01
t = 0.1
Nt = int(t/dt)+1

folder = "Diffusion"
filename = "time"
Timeseries = []

color = ['b', 'g', 'r','c', 'm', 'k']
cc = 0
for i in range(0,Nt,5):
   #Calculate the analytic solution for
    x = np.arange(-1,1,0.015)
    t = i*dt
    if t == 0:
        A = IC(x,t)
    else:
        A = G(x,t)
    ana, = plt.plot(x,A,'ro', label = 'Analytical Solution')

    try:
        Data = np.loadtxt(folder+"/"+filename+"_%d.txt"%(i))
        sim, = plt.plot(Data[:,0], Data[:,2],color='blue', linewidth = 2.0, label ='Simulation' )
    # Plot your Data 
        #print 'Sim'
        #print np.var(Data[:,2])
    #Timeseries.append(Data)
    except: pass
    cc += 1
#Timeseries = np.array(Timeseries)
#ana = mpatches.Patch(color='red', label='')
#sim = mpatches.Patch(color='red', label='Simulation')

plt.legend([ana, sim], [ 'Analytical Solution','Simulation' ])
plt.ylabel('Concentration')
plt.xlabel('x-coordinate')
plt.show()


