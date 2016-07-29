import numpy as np 
import matplotlib.pyplot as plt

#Particle with postion concentration and velocity
class particle:
   def __init__(self,x,v,c):
       self.x = x
       self.v = v
       self.c = c 
       self.NN = []

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
folder = "Error"
filename = "h01t001"

#Discretization in time
dt = 0.1
t = 1.0
ntMax = int(t/dt)+1


#Discretization in spce
h = 0.1
L = 2.0

#Kernel 
eps = h*4.0
Vq = h

#Cutoff
rc = eps*2.0

#Parameters:
PiAD = 1.0
PiA1 = 1.0
PiA2 = 1.0
PiA3 = 1.0
PiA4 = 1.0

PiPD = 1.0
PiP1 = 1.0
PiP2 = 1.0
PiP3 = 1.0
PiP4 = 1.0

alpha = 2.0
beta  = 2.0


#Create a list of prticles on regular Grid betwee  x = -1 and 1 

# Write the initial timestep
with file(folder+"/timestep_%d.txt"%(0), "w") as outputfile:
   N = int(L/h)
   ParticleList = []
   for i in range(N):
      x = -1.0 + i*h
      v = .0
       #c = np.array([0.0,0.0])
      c = np.array([1.0,0.0])*np.exp(-x**2/0.05)+0.1
      p = particle(x,v,c)
      ParticleList.append(p)
 
   for p in ParticleList:
      outputfile.write("%0.8f %0.8f %0.8f %0.8f \n"%(p.x,p.v,p.c[0],p.c[1])) 


# Integration 
for i in range(1,ntMax):
    #Cell List 
    for p in ParticleList:
        p.NN = []
        for q in ParticleList:
           dx = q.x-p.x
           if np.abs(dx) > 0 and np.abs(dx) < rc:
              p.NN.append(q)

           if  np.abs(dx-L) < rc+h/2:
              x = q.x-L
              Ghost = particle(x,q.v,q.c)
              p.NN.append(Ghost)

           if np.abs(dx+L) < rc:
              x = q.x + L            
              Ghost = particle(x,q.v,q.c)
              p.NN.append(Ghost)
        #print p.NN
              
    #Calculate the new particle postions
    flag = True
    for p in ParticleList:

        #Change in concentration

        dc = np.array([.0, .0])

        #Transport
        for q in p.NN:
            dx = q.x - p.x
            #Diffusion
            dc[0] += PiAD/eps**2*Vq*(q.c[0]-p.c[0])*O2_1Dkernel((dx)/eps)/eps 
            dc[1] += PiPD/eps**2*Vq*(q.c[1]-p.c[1])*O2_1Dkernel((dx)/eps)/eps 
            #if flag:
            #   print O2_1Dkernel((dx)/eps)/eps
            #   print dx

        #if flag:
        #   print '-------------'
        #flag = False 
            #Advection 
            dc[0] -= p.c[0]/eps*Vq*(q.v + p.v)*O1_1Dkernel((dx)/eps)/eps 
            dc[1] -= p.c[1]/eps*Vq*(q.v + p.v)*O1_1Dkernel((dx)/eps)/eps

        #Reactions
        #dc[0] -= ((PiA1+PiA3)*p.c[0] + PiA4 * p.c[1]**alpha * p.c[0] - PiA2)
        #dc[1] -= ((PiP1+PiP3)*p.c[1] + PiP4 * p.c[0]**beta  * p.c[1] - PiP2)

        #Integrate the concentrations
        p.c[0] += dc[0]*dt
        p.c[1] += dc[1]*dt

        #Equation of motion:
        p.x += p.v*dt

        #Apply Periodic boundary conditions 
        if p.x > L/2:
            p.x = p.x - L
        if p.x < -L/2:
            p.x = p.x + L
    
    #Create a new particle list and interpolate the quantities from the old particles 
    #Write the postions, the velocities and the concentrations in an output file
            
    with file(folder+"/"+filename+"_%d.txt"%(i), "w") as outputfile:
       NewParticleList = []
       for k in range(N):
 
          x = -1.0 + k*h
          v = .0
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


        


    
