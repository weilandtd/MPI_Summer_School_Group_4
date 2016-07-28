#using PyCall
using Plots
using Sundials

#Parameters:

A = 1.0
P = 1.0

kAP = 1.0
kPA = 1.0
konA = 1.0
koffA = 1.0
konP = 1.0
koffP = 1.0
NA = 1,0
V = 1.0
Amem = 1.0
alpha = 1.0
beta = 1.0

dt = 0.2
t = 0.0

max = 1000
t_max = max*dt

abstract p

function f(t, y, ydot,p)
    ydot[1] = y[1]*p.kA1 - p.kA2*y[2]^p.alpha*y[1] - p.kA3 
    ydot[2] = y[2]*p.kP1 - p.kP2*y[1]^p.beta *y[2] - p.kP3
end

t = [0.0, 4 * logspace(-1., 7., 9)]
res = Sundials.cvode(f, [1.0, 0.0, 0.0], t)



plot(t_,AP[:,1], label = "aPAR", color="blue" );
plot!(t_,AP[:,2], label = "pPAR"  , color="red", xlabel="time", ylabel="c(t)" )

#plt.xlim([0,dt*max])

savefig("ODE_SystemPAR.svg")
