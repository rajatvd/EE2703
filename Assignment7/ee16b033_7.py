
import sys
import argparse
from pylab import *

"""
Get command line arguments and assign them to variables automatically.
"""

ARGS = {}

def addargs(args,parser):
    for arg in args:
        parser.add_argument('-{}'.format(arg),default=args[arg],type=type(args[arg]))
        
ap = argparse.ArgumentParser()
addargs(ARGS,ap)
args = ap.parse_args()

for arg in ARGS:
    exec("{} = args.{}".format(arg,arg))
    

import scipy.signal as sp

def F_s(freq=1.5,decay=0.5):
    """Transfer function of the given system"""
    n = poly1d([1,decay])
    d = n*n+freq**2
    return n,d

def secondOrderH(wn=1.5,zeta=0,gain=1/2.25):
    """General second order all pole transfer function."""
    n = poly1d([wn**2*gain])
    d = poly1d([1,2*wn*zeta,wn**2])
    return n,d

def solveProblem(decay):
    """Find the response to the given system to a decaying cosine."""
    inN, inD = F_s(decay=decay)
    HN, HD = secondOrderH()

    outN,outD = inN*HN, inD*HD

    out_s = sp.lti(outN.coeffs, outD.coeffs)

    t = linspace(0,50,1000)
    return sp.impulse(out_s,None,t)

figure()
plot(*solveProblem(0.5))
title("System response with decay of $0.5$")
xlabel("$t$")
ylabel("$x$")
show()

figure()
plot(*solveProblem(0.05))
title("System response with decay of $0.05$")
xlabel("$t$")
ylabel("$x$")
show()

def input_f(t,decay=0.5,freq=1.5):
    """Exponentially decaying cosine function."""
    u_t = 1*(t>0)
    return cos(freq*t)*exp(-decay*t) * u_t
    

# get transfer function
system = secondOrderH()

# time range
t = linspace(0,70,1000)

# list of outputs
outs = []

# list of frequencies to iterate over
fs = linspace(1.4,1.6,5)


for freq in fs:
    # solve
    t,y,svec = sp.lsim(system,input_f(t,decay=0.05,freq=freq),t)
    
    # store
    outs.append(y)

# plot stored outputs
figure()
plot(t,array(outs).T,linewidth=0.7)
title("System responses with variation of input frequency")
xlabel("$t$")
ylabel("$x$")
legend(["Freq = ${:.2f}$".format(f) for f in fs], loc=2)
show()

w,S,phi=sp.lti(*system).bode()

semilogx()
plot(w,S)
title("Magnitude plot")
xlabel("Frequency in rad/s (log)")
ylabel("Magnitude in dB")
grid()
show()

semilogx()
plot(w,phi)
title("Phase plot")
xlabel("Frequency in rad/s (log)")
ylabel("Phase in degrees")
grid()
show()

X_s = sp.lti([1,0,2],[1,0,3,0])

Y_s = sp.lti([2],[1,0,3,0])

t = linspace(0,20,1e3)

t, x = sp.impulse(X_s,None,t)
t, y = sp.impulse(Y_s,None,t)

figure()
title("Responses of coupled system")
ylabel("Displacement")
xlabel("$t$")
plot(t,x)
plot(t,y)
legend(["$x(t)$", "$y(t)$"])
show()

# Find the transfer function of the given circuit
R = 100
L = 1e-6
C = 1e-6

wn = 1/sqrt(L*C) # natural frequency
Q = 1/R * sqrt(L/C) # quality factor
zeta = 1/(2*Q) # damping constant

# transfer function
n,d = secondOrderH(gain=1,wn=wn,zeta=zeta)

# make system
H = sp.lti(n,d)

# get bode plots
w,S,phi=H.bode()

semilogx()
plot(w,S)
title("Magnitude plot")
xlabel("Frequency in rad (log)")
ylabel("Magnitude in dB")
grid()
show()

semilogx()
plot(w,phi)
title("Phase plot")
xlabel("Frequency in rad (log)")
ylabel("Phase in degrees")
grid()
show()

def input_2(t,w1=1e3,w2=1e6):
    """Two cosines of different frequencies"""
    u_t = 1*(t>0)
    return (cos(w1*t)-cos(w2*t)) * u_t

# early response
t1=linspace(0,30e-6,1e3)
t1,y1,svec = sp.lsim(H,input_2(t1),t1)

# steady state response
t2=linspace(0,10e-3,1e3)
t2,y2,svec = sp.lsim(H,input_2(t2),t2)

figure()
title(r"Response for 30 $\mu s$")
xlabel("$t$ (sec)")
ylabel("$v_0(t)$")
plot(t1,y1)
show()

figure()
title(r"Response for 10 msec")
xlabel("$t$ (sec)")
ylabel("$v_0(t)$")
plot(t2,y2)
show()