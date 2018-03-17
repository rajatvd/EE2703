
import sys
import argparse

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
    n = poly1d([1,decay])
    d = n*n+freq**2
    return n,d

def secondOrderH(wn=1.5,zeta=0,gain=1/2.25):
    n = poly1d([wn**2*gain])
    d = poly1d([1,2*wn*zeta,wn**2])
    return n,d

def solveProblem(decay):
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
    out = zeros(t.shape)
    ii = where(t>0)
    out[ii] = cos(freq*t[ii])*exp(-decay*t[ii])
    return out

system = secondOrderH()
t = linspace(0,50,1000)
outs = []
fs = linspace(1.4,1.6,5)
for freq in fs:
    t,y,svec = sp.lsim(system,input_f(t,decay=0.05,freq=freq),t)
    outs.append(y)

figure()
plot(t,array(outs).T,linewidth=0.7)

title("System responses with variation of input frequency")
xlabel("$t$")
ylabel("$x$")

legend(["Freq = ${:.2f}$".format(f) for f in fs], loc=2)

show()

figure()
plot(t,array(outs).T,linewidth=0.7)

title("System responses with variation of input frequency")
xlabel("$t$")
ylabel("$x$")
legend(["Freq = ${:.2f}$".format(f) for f in fs], loc=2)
show()

R = 100
L = 1e-6
C = 1e-6

wn = 1/sqrt(L*C)
Q = 1/R * sqrt(L/C)
zeta = 1/(2*Q)

n,d = secondOrderH(gain=1,wn=wn,zeta=zeta)