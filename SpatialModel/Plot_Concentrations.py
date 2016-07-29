import numpy as np 
import matplotlib.pyplot as plt

dt = 0.01
t = 1.0
ntMax = int(t/dt)+1

folder = "Diffusion"
filename = "time"
Timeseries = []
for i in range(Nt,0,10):
    Data = np.loadtxt(folder+"/"+filename+"_%d.txt"%(i))
    plt.plot(Data[:,0], Data[:,2] )
    #Timeseries.append(Data)
    

for i

#Timeseries = np.array(Timeseries)
plt.show()


