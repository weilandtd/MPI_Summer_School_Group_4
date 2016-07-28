import numpy as np 
from matplotlib import plot as plt


#Particle with postion concentration and velocity
class particle:
   def __init__(x,v,c):
       particle.x = x
       particle.v = v
       particle.c = c 
       NN = []

#1D Kernel
def 1Dkernel(x):
    return 
   
#Parameters
class parameters:
    def __init__(x,v,c):


#Create a list of prticles on regular Grid betwee  x = -1 and 1 
h = 0.2
N = np.floor(2/h)
ParticleList = []
for i in range(N):
    x = -1.0 + i*h
    v = 0.1
    c = np.array([1.0,0.0])
    ParticleList.append(particle(x,v,c))
    

#Simulate 
dt = 0.1
ntMax = 1000
for i in range(ntMax):
    for particle in ParticleList:
        #Equation of motion:
        particle.x += particle.v*dt
        for q in NN:
             += 
             +=
 
        # Periodic boundary conditions 
        if particle.x > 1.0:
            particle.x = particle.x - 2.0
        if particle.x < -1.0:
            particle.x = particle.x + 2.0

        with file = open("timestep_%d.txt"%(i), "w"):
            file.write("%0.8f %0.8f %0.8f"%(particle.x,particle.v,particle.c[1],))

