
import sys
import argparse
from pylab import *
from IPython.display import *

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
    

from numpy.fft import *
rcParams['figure.figsize'] = 18,6
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
x=rand(10)
X= fft(x)
y= ifft(X)
print("Original and reconstructed signal values:")
print(c_[x,y])
print("Max absolute error in reconstruction: {}".format(abs(x-y).max()))

x=linspace(0,2*pi,128)
y=sin(5*x)
Y= fft(y)
figure()
subplot(1,2,1)
plot(abs(Y),lw=2)
ylabel(r"$|Y|$",size=16)
xlabel(r"$k$",size=16)
title(r"Magnitude of FFT of $\sin(5t)$")
grid(True)
subplot(1,2,2)
plot(unwrap(angle(Y)),lw=2)
ylabel(r"Phase of $Y$",size=16)
title(r"Phase of FFT of $\sin(5t)$")
xlabel(r"$k$",size=16)
grid(True)
show()

def plotFFT(func,t_range=(0,2*pi),points=128,tol=1e-5,
            func_name=None,unwrap_=True,wlim=(-10,10),scatter_size=40,
           iff=False):
    """Plot the FFT of the given continuous function.
    
    func : the continuous function
    t_range : the time range over which to sample the function,
        exclusive of the last value
    points : number of samples
    tol : tolerance for setting phase to 0 when magnitude is low
    func_name : name of the function
    unwrap : whether to unwrap phase
    wlim : range of frequencies for the plots, give None for all frequencies
    scatter_size : size of scatter plot points
    iff: whether to do an ifftshift on the time range
    
    Returns:
    numpy array containing the FFT, after being shifted and normalized.
    """
    
    # default name for function
    if func_name == None:
        func_name = func.__name__
    
    # time points to sample
    t = linspace(*t_range,points+1)[:-1]
    T = t_range[1]-t_range[0]
    samplingfreq = points/T
    
    if iff:
        t = ifftshift(t)
    
    # corresponding frequencies of the sampled signal
    w = linspace(-pi,pi,points+1)[:-1]
    w = w*samplingfreq
    
    # find fft
    y = func(t)
    Y =  fftshift( fft(y))/points
    
    # get phase
    ph = angle(Y)
    if unwrap_:
        ph = unwrap(ph)
    
    # get mag
    mag = abs(Y)
    
    # clean up phase where mag is sufficiently close to 0
    ph[where(mag<tol)]=0
    
    # plot 
    fig,axes = subplots(1,2)
    ax1,ax2 = axes
    
    # magnitude
    ax1.set_title("Magnitude of DFT of {}".format(func_name))
    ax1.set_xlabel("Frequency in rad/s")
    ax1.set_ylabel("Magnitude")
    ax1.plot(w,mag,color='red')#,s=scatter_size)
    ax1.set_xlim(wlim)
    ax1.grid()
    
    # phase
    ax2.set_title("Phase of DFT of {}".format(func_name))
    ax2.set_xlabel("Frequency in rad/s")
    ax2.set_ylabel("Phase in rad")
    ax2.scatter(w,ph,color='green',s=scatter_size)
    ax2.set_xlim(wlim)
    ax2.grid()
    
    show()
    return w,Y

def sin5(t):
    return sin(5*t)

f = plotFFT(sin5,unwrap_=False, func_name='$sin(5t)$')

t=linspace(0,2*pi,129);t=t[:-1]
y=(1+0.1*cos(t))*cos(10*t)
Y= fftshift( fft(y))/128.0
w=linspace(-64,63,128)
figure()
subplot(1,2,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r"$|Y|$",size=16)
xlabel(r"$\omega$",size=16)
title(r"Magnitude of FFT of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$")
grid(True)
subplot(1,2,2)
title(r"Phase of FFT of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$")
plot(w,angle(Y),'ro',lw=2)
xlim([-15,15])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

def amMod(t):
    return (1+0.1*cos(t))*cos(10*t)

f = plotFFT(amMod,unwrap_=False,wlim=(-15,15),t_range=(-4*pi,4*pi),points=512,func_name=r"$(1+0.1\cos(t))\cos(10t)$")

def sin3(t):
    return sin(t)**3

f=plotFFT(sin3,unwrap_=False,t_range=(-4*pi,4*pi),points=512,wlim=(-5,5),func_name=r"$\sin^3(t)$")

def cos3(t):
    return cos(t)**3

f=plotFFT(cos3,unwrap_=False,t_range=(-4*pi,4*pi),points=512,wlim=(-5,5),func_name=r"$\cos^3(t)$")

def fm(t):
    return cos(20*t + 5*cos(t))

f=plotFFT(fm,wlim=(-40,40),unwrap_=False,t_range=(-16*pi,16*pi),points=2**10,func_name=r"$ \cos(20t +5 \cos(t))$")

def gauss(t):
    return exp(-t**2/2)

def expected_gaussFT(w):
    return 1/sqrt(2*pi) * exp(-w**2/2)

def estimateCTFT(func, tol=1e-6,time_samples=128, true_func=None, 
                 func_name=None, wlim=None, scatter_size=40):
    """Estimate the continuous time Fourier Transform of the given function
    by finding the DFT of a sampled window of the function. The magnitude and
    phase of the estimate are also plotted.
    
    The window size and sample number are doubled until the consecutive 
    total absolute error between two estimates is less than the given tolerance.
    
    time_samples : the initial number of samples to start with
    true_func : A function which is the analytical CTFT of the given function.
                Used to compare the estimate results with the true results.
    
    Returns frequencies and the CTFT estimate at those frequencies.
    """
    
    if func_name == None:
        func_name = func.__name__
    
    
    T = 8*pi
    N = time_samples
    Xold = 0
    error = tol+1
    iters=0
    
    while error>tol:
        
        delta_t = T/N # time resolution
        delta_w = 2*pi/T # frequency resolution

        W = N*delta_w # total frequency window size

        t = linspace(-T/2,T/2,N+1)[:-1] # time points
        w = linspace(-W/2,W/2,N+1)[:-1] # freq points

        x = func(t)

        # find DFT and normalize
        # note that ifftshift is used to prevent artifacts in the
        # phase of the result due to time domain shifting
        X = delta_t/(2*pi) * fftshift(fft(ifftshift(x)))
        
        error = sum(abs(X[::2]-Xold))
        
        Xold = X
        N *= 2 # number of samples
        T *= 2 # total time window size
        iters+=1
        
        
    print("Estimated error after {} iterations: {}".format(iters, error))
    print("Time range : ({:.4f}, {:.4f})".format(-T/2,T/2))
    print("Time resolution : {:.4f}".format(delta_t))
    print("Frequency resolution : {:.4f}".format(delta_w))
        
    if true_func != None:
        true_error = sum(abs(X-true_func(w)))
        print("True error: {}".format(true_error))
    
    mag = abs(X)
    ph = angle(X)
    ph[where(mag<tol)]=0
    
    # plot estimate
    fig,axes = subplots(1,2)
    ax1,ax2 = axes
    
    # magnitude
    ax1.set_title("Magnitude of CFT estimate of {}".format(func_name))
    ax1.set_xlabel("Frequency in rad/s")
    ax1.set_ylabel("Magnitude")
    ax1.plot(w,mag,color='red')#,s=scatter_size)
    ax1.set_xlim(wlim)
    ax1.grid()
    
    # phase
    ax2.set_title("Phase of CFT estimate of {}".format(func_name))
    ax2.set_xlabel("Frequency in rad/s")
    ax2.set_ylabel("Phase in rad")
    ax2.scatter(w,ph,color='green',s=scatter_size)
    ax2.set_xlim(wlim)
    ax2.grid()
    
    show()
    
    if true_func != None:
        
        X_ = true_func(w)
        
        mag = abs(X_)
        ph = angle(X_)
        ph[where(mag<tol)]=0
        
        # plot true
        fig,axes = subplots(1,2)
        ax1,ax2 = axes

        # magnitude
        ax1.set_title("Magnitude of true CFT of {}".format(func_name))
        ax1.set_xlabel("Frequency in rad/s")
        ax1.set_ylabel("Magnitude")
        ax1.plot(w,mag,color='red')#,s=scatter_size)
        ax1.set_xlim(wlim)
        ax1.grid()

        # phase
        ax2.set_title("Phase of true CFT {}".format(func_name))
        ax2.set_xlabel("Frequency in rad/s")
        ax2.set_ylabel("Phase in rad")
        ax2.plot(w,ph,color='green')
        ax2.set_xlim(wlim)
        ax2.grid()
        
        show()
    return w,X

w_,X_ = estimateCTFT(gauss,true_func=expected_gaussFT,wlim=(-10,10),func_name=r"$e^{\frac{-t^2}{2}}$")