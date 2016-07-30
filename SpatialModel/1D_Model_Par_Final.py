import numpy as np 
import time  
import matplotlib.pyplot as plt

#Particle with postion concentration and velocity
class particle:
   def __init__(self,x,v,c, l = 0):
       self.x = x
       self.v = v
       self.c = c 
       self.NN = []
       self.l = l
       self.dc = np.array([.0 ,.0])

#1D Kernel kern
def O2_1Dkernel(x):
    return 1.0/(np.sqrt(np.pi))*np.exp(-x**2)

def O1_1Dkernel(x):
    return x/(np.sqrt(np.pi))*np.exp(-x**2)

#Interpolation function 
def Interpolation(s):
    if s < 1.0:
        return 1.0-0.5*(5.0*s**2-3.0*s**3)
    elif s <= 2.0:
       return 0.5*(2.0 -s )**2 * (1-s)
    else: return .0

#OutputFile
folder = "PAR"
filename = "time"

#Discretioation in time
dt = 0.02
# 550s / k_offP
#t = 76495.0
t = 10.0
ntMax = int(t/dt)+1

#Discretization in spce
h = 0.02
L = 2.0

#Kernel 
eps = 0.2 #h*15.0
Vq = h*1.0

#Calculate a lookup for ther Kernel
h_ = np.arange(0,500,1.0)*h 
print h_

eta_2_eps =  O2_1Dkernel(h_/eps)/eps
N_ = 500
h_ = np.arange(-N_,N_+1,1.0)*h 
eta_1_eps =  O1_1Dkernel(h_/eps)/eps

print eta_2_eps

#Cutoff
rc = 3.0*eps

#In real space:
V = 2.5*10**4
omega = 4.4*10**3
N_A = 2.4*10**5
N_P = 9.8*10**4
D_A = 0.28
D_P = 0.15
k_offA = 3.25*10**-3
k_offP = 7.19*10**-3
k_onA  = 6.29*10**-3
k_onP  = 7.682*10**-2
k_AP = 0.004
k_PA = 0.01

alpha = 2.0
beta  = 2.0

#Dimensionless Parameters:
#Created on Sat Jul 30 14:14:18 2016
#@author: Mahdiye
# Dimensionless
#omega = 1.0
#ro = gamma = 1/omega
#tau = 1/k_offP

PiAD = D_A/(k_offP*omega)
PiA1 = k_offA/k_offP
PiA2 = k_onA*omega*N_A/(k_offP*V)
PiA3 = k_onA*omega/(k_offP*V)
PiA4 = k_AP/(k_offP*omega**alpha)

PiPD = D_P/(k_offP*omega)
PiP1 = 1
PiP2 = k_onP*omega*N_P/(k_offP*V)
PiP3 = k_onP*omega/(k_offP*V)
PiP4 = k_PA/(k_offP*omega**beta)



#Create a list of prticles on regular Grid betwee  x = -1 and 1 
D = 1.0 
t0 = 0.001

# Write the initial timestep
with file(folder+"/"+filename+"_%d.txt"%(0), "w") as outputfile:
   N = int(L/h)
   ParticleList = []
   for i in range(0,N):
      x = -1.0 + i*h
      v = 0.0
       #c = np.array([0.0,0.0])
      c = np.array([0.0,10000.0])
      p = particle(x,v,c)
      ParticleList.append(p)
 
   for p in ParticleList:
      outputfile.write("%0.8f %0.8f %0.8f %0.8f \n"%(p.x,p.v,p.c[0],p.c[1])) 

N_index = int((rc/h))

# Integration 
for i in range(1,ntMax):
    #Cell List 
    for n in range(0,N):
        p = ParticleList[n]
        p.NN = []
        # Neigbours
        P_index = n #np.int_((p.x+L/2.0)/h)
        #print (p.x+L/2.0)/h
        #print P_indexxs
        #print '--------------'
        
        for k in range(P_index-N_index, P_index+N_index+1,1):
           if k == P_index:
              pass
           elif k >= 0 and k < N:
              index = k
              x = ParticleList[index].x
              v = ParticleList[index].v
              c = ParticleList[index].c
              l = k - P_index
              p.NN.append(particle(x,v,c,l))
           elif k >= N:
              index = k - N 
              iq = 1.0          
              while index >= N :
                 index = index-N
                 iq += 1.0
              x = ParticleList[index].x + L*iq
              v = ParticleList[index].v
              c = ParticleList[index].c
              l = k - P_index
              p.NN.append(particle(x,v,c,l))
           elif k < 0:
              index = k + N 
              iq = 1.0
              while index < 0:
                 index = index + N
                 iq += 1.0

              x = ParticleList[index].x - L*iq
              v = ParticleList[index].v
              c = ParticleList[index].c
              l = k - P_index
              p.NN.append(particle(x,v,c,l))


        #for q in ParticleList:
        #   dx = q.x-p.x
        #   if np.abs(dx) > 0 and np.abs(dx) < rc:
        #      p.NN.append(q)
        #   
        #   if  np.abs(dx) > 0 and  np.abs(dx-L) < rc :
        #      x = q.x-L
        #      Ghost = particle(x,q.v,q.c)
        #      p.NN.append(Ghost)
        #
        #   if  np.abs(dx) > 0 and np.abs(dx+L) < rc :
        #      x = q.x + L            
        #      Ghost = particle(x,q.v,q.c)
        #      p.NN.append(Ghost)
              
    #Calculate the new particle quantity
    for p in ParticleList:
       #if p.x == 0.0:
       #   flag = True

       test = 0.0
        #Change in concentration
       dc = np.array([.0, .0])
       blub = 0.0
        #Transport
       for q in p.NN:
          dx = q.x - p.x

          #Diffusion
          #print q.l
          l = abs(q.l)
          #print l
          dc[0] += PiAD/(eps**2)*Vq*(q.c[0]-p.c[0])*eta_2_eps[l]
          dc[1] += PiPD/(eps**2)*Vq*(q.c[1]-p.c[1])*eta_2_eps[l]
          
          blub += dx*eta_2_eps[l]
          #Advection 
          dc[0] -= p.c[0]/eps*Vq*(q.v + p.v)*eta_1_eps[q.l+N_/2]
          dc[1] -= p.c[1]/eps*Vq*(q.v + p.v)*eta_1_eps[q.l+N_/2]

        #Reactions
       dc[0] -= ((PiA1+PiA3)*p.c[0] + PiA4 * p.c[1]**alpha * p.c[0] - PiA2)
       dc[1] -= ((PiP1+PiP3)*p.c[1] + PiP4 * p.c[0]**beta  * p.c[1] - PiP2)
       
       p.dc = dc

       #print blub
       #print '--------------'

    #Integrate the concentrations
    for p in ParticleList:
       p.c[0] += p.dc[0]*dt
       p.c[1] += p.dc[1]*dt

        #Equation of motion:
       p.x += p.v*dt

        #Apply Periodic boundary conditions 
       if p.x > L/2.0:
          p.x = p.x - L
       if p.x < -L/2.0:
          p.x = p.x + L
       
    
    #Create a new particle list and interpolate the quantities from the old particles 
    #Write the postions, the velocities and the concentrations in an output file
            
    with file(folder+"/"+filename+"_%d.txt"%(i), "w") as outputfile:
       NewParticleList = []
       for k in range(0,N):
 
          x = -1.0 + k*h
          v = 0.0
          c = np.array([.0 , .0])
          for p in ParticleList:
             dx = x - p.x

             if dx > L/2:
                dx = dx - L
             if dx < -L/2:
                dx = dx + L
                  
             s = abs(dx/h)
             c += p.c*Interpolation(s)
             
    
          #Write out the data from the grid
          NewParticleList.append(particle(x,v,c))
          outputfile.write("%0.8f %0.8f %0.8f %0.8f \n"%(x,v,c[0],c[1])) 
       outputfile.write("\n") 

  
       #for p in ParticleList:    
          #outputfile.write("%0.8f %0.8f %0.8f %0.8f \n"%(p.x,p.v,p.c[0],p.c[1])) 


    del ParticleList
    #Assing new Particles 
    ParticleList = NewParticleList


        


    
