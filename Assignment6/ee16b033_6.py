import sys
import argparse

# Get command line arguments
ARGS = {"n":100,
       "M":5,
       "nk":500,
       "u0":5,
       "p":0.25,
       "Msig":1}

def addargs(args,parser):
    for arg in args:
        parser.add_argument('-{}'.format(arg),default=args[arg],type=type(args[arg]))
        
ap = argparse.ArgumentParser()
addargs(ARGS,ap)
args = ap.parse_args()

for arg in ARGS:
    exec("{} = args.{}".format(arg,arg))
    

    

from pylab import *
# Increase figure and font size
rcParams['figure.figsize'] = 12,9
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

def simulateTubelight(n,M,nk,u0,p,Msig):
    """
    Simulate a tubelight and return the electron positions and velocities,
    and positions of photon emissions.
    
    n: integer length of tubelight
    M: average number of electrons generated per timestep
    nk: total number of timesteps to simulate
    u0: threshold voltage for ionization
    p: probability of ionization given an electron is faster than the threshold
    Msig: stddev of number of electrons generated per timestep
    
    """

    xx = zeros(n*M)
    u = zeros(n*M)
    dx = zeros(n*M)

    I = []
    X = []
    V = []

    for k in range(nk):

        # add new electrons
        m=int(randn()*Msig+M)
        jj = where(xx==0)
        xx[jj[0][:m]]=1

        # find electron indices
        ii = where(xx>0)

        # add to history lists
        X.extend(xx[ii].tolist())
        V.extend(u[ii].tolist())

        # update positions and speed
        dx[ii] = u[ii]+0.5
        xx[ii]+=dx[ii]
        u[ii]+=1

        # anode check
        kk = where(xx>=n)
        xx[kk]=0
        u[kk]=0

        # ionization check
        kk = where(u>u0)[0]
        ll=where(rand(len(kk))<=p);
        kl=kk[ll];

        # ionize
        dt = rand(len(kl))
        #xx[kl]=xx[kl]-dx[kl]+((u[kl]-1)*dt+0.5*dt*dt)
        xx[kl]=xx[kl]-dx[kl]*dt

        u[kl]=0

        # add emissions
        I.extend(xx[kl].tolist())
        
    return X,V,I



def plotGraphs(X,V,I):
    """
    Plot histograms for X and I, and a phase space using X and V.
    Returns the emission intensities and locations of histogram bins.
    """
    
    # electron density
    figure()
    hist(X,bins=n,cumulative=False)
    title("Electron density")
    xlabel("$x$")
    ylabel("Number of electrons")
    show()

    # emission instensity
    figure()
    ints,bins,trash = hist(I,bins=n)
    title("Emission Intensity")
    xlabel("$x$")
    ylabel("I")
    show()

    # electron phase space
    figure()
    scatter(X,V,marker='x')
    title("Electron Phase Space")
    xlabel("$x$")
    ylabel("$v$")
    show()
    
    return ints,bins


X,V,I = simulateTubelight(n,M,nk,u0,p,Msig)
ints, bins = plotGraphs(X,V,I)

# Tabulate emission counts
xpos=0.5*(bins[0:-1]+bins[1:])
from tabulate import *
print("Intensity Data")
print(tabulate(stack((xpos,ints)).T,["xpos","count"]))
