
import sys
import argparse
from pylab import *
from IPython.display import *
rcParams['figure.figsize'] = 18,6
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

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
    

t=linspace(-pi,pi,65);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
y=sin(sqrt(2)*t)
y[0]=0 # the sample corresponding to -tmax should be set zeroo
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

t1=linspace(-pi,pi,65);t1=t1[:-1]
t2=linspace(-3*pi,-pi,65);t2=t2[:-1]
t3=linspace(pi,3*pi,65);t3=t3[:-1]
# y=sin(sqrt(2)*t)
figure(2)
plot(t1,sin(sqrt(2)*t1),"b",lw=2)
plot(t2,sin(sqrt(2)*t2),"r",lw=2)
plot(t3,sin(sqrt(2)*t3),"r",lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)$")
grid(True)
show()

t1=linspace(-pi,pi,65);t1=t1[:-1]
t2=linspace(-3*pi,-pi,65);t2=t2[:-1]
t3=linspace(pi,3*pi,65);t3=t3[:-1]
y=sin(sqrt(2)*t1)
figure(3)
plot(t1,y,"bo",lw=2)
plot(t2,y,"ro",lw=2)
plot(t3,y,'ro',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)$ with $t$ wrapping every $2\pi$ ")
grid(True)
show()

t=linspace(-pi,pi,65);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
y=t
y[0]=0 # the sample corresponding to -tmax should be set zeroo
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
figure()
semilogx(abs(w),20*log10(abs(Y)),lw=2)
xlim([1,10])
ylim([-20,0])
xticks([1,2,5,10],["1","2","5","10"],size=16)
ylabel(r"$|Y|$ (dB)",size=16)
title(r"Spectrum of a digital ramp")
xlabel(r"$\omega$",size=16)
grid(True)
show()

def hamming(n):
    n = array(n)
    N = n.shape[0]
    return 0.54+0.46*cos(2*pi*n/(N-1))

n = arange(64)
plot(n,fftshift(hamming(n)))
scatter(n,fftshift(hamming(n)))
grid()
xlim(0,63)
ylabel("$y$")
xlabel("$n$")
title(r"Hamming window $0.54+0.46\cos(\frac{2\pi n}{N-1})$ for $N = 64$")
show()


t1=linspace(-pi,pi,65);t1=t1[:-1]

t2=linspace(-3*pi,-pi,65);t2=t2[:-1]
t3=linspace(pi,3*pi,65);t3=t3[:-1]
n=arange(64)
wnd=fftshift(0.54+0.46*cos(2*pi*n/63))
y=sin(sqrt(2)*t1)*wnd
figure(3)
plot(t1,y,'bo',lw=2)
plot(t2,y,'ro',lw=2)
plot(t3,y,'ro',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)\times w(t)$ with $t$ wrapping every $2\pi$ ")
grid(True)
show()

t=linspace(-pi,pi,65);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
n=arange(64)
wnd=fftshift(0.54+0.46*cos(2*pi*n/63))
y=sin(sqrt(2)*t)*wnd
y[0]=0 # the sample corresponding to -tmax should be set zeroo
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-8,8])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-8,8])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

t=linspace(-4*pi,4*pi,257);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
n=arange(256)
wnd=fftshift(0.54+0.46*cos(2*pi*n/256))
y=sin(sqrt(2)*t)
# y=sin(1.25*t)
y=y*wnd
y[0]=0 # the sample corresponding to -tmax should be set zeroo
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/256.0
w=linspace(-pi*fmax,pi*fmax,257);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),"b",w,abs(Y),"bo",lw=2)
xlim([-4,4])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),"ro",lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y$",size=16)
12
xlabel(r"$\omega$",size=16)
grid(True)
show()

def plotFFT(func,t_range=(0,2*pi),points=128,tol=1e-5,
            func_name=None,unwrap_=True,wlim=(-10,10),scatter_size=40,
           iff=False, plot=True, window=False):
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
    
    y = func(t)
    if iff:
        y[0]=0
        y = fftshift(y)
    
    # corresponding frequencies of the sampled signal
    w = linspace(-pi,pi,points+1)[:-1]
    w = w*samplingfreq
    
    # find fft
    if window:
        wnd = fftshift(hamming(arange(points)))
        y = y*wnd
    Y =  fftshift( fft(y))/points
    
    if not plot:return w,Y
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
    ax1.plot(w,mag,color='red')
    ax1.scatter(w,mag,color='red',s=scatter_size)
    ax1.set_xlim(wlim)
    ax1.grid()
    
    # phase
    ax2.set_title("Phase of DFT of {}".format(func_name))
    ax2.set_xlabel("Frequency in rad/s")
    ax2.set_ylabel("Phase in rad")
    ax2.plot(w,ph,color='green')
    ax2.scatter(w,ph,color='green',s=scatter_size)
    ax2.set_xlim(wlim)
    ax2.grid()
    
    show()
    return w,Y

w,Y = plotFFT(lambda x: cos(0.86*x)**3, points = 2048, t_range=(-4*pi,4*pi), iff=True,
              wlim=(-4,4),tol=1e-2, window=False,func_name=r"$\cos^3(0.86t)$",unwrap_=False)

w,Y = plotFFT(lambda x: cos(0.86*(x))**3, points = 2048, t_range=(-4*pi,4*pi), 
              wlim=(-4,4),tol=1e-2, window=True,func_name=r"$\cos^3(0.86t)$",unwrap_=False)

mag = abs(Y)
l = 5 # must be odd
index = argmax(mag)
ind = index-int((l-1)/2)+ arange(l)
fs = w[ind]
points = mag[index-int((l-1)/2):index+int((l+1)/2)]
freq = np.sum(points*fs)/np.sum(points)
print("Estimated first peak frequency: {:.4f}".format(abs(freq)))


def estimateWD1(vec,l=6,window=True):
    """Estimate the value of omega and delta assuming vec contains
    128 samples of cos(omega*t + delta) in (-pi,pi). Uses the magnitude
    and phase spectra to estimate omega and delta respectively."""
    
    N = 128
    delta_t = 2*pi/N
    w_max = pi/delta_t
    delta_w = 2*w_max/N
    
    w = linspace(-w_max,w_max,N+1)[:-1]
    
    if window:
        vec_ = vec*fftshift(hamming(arange(N)))
    else:
        vec_ = vec
        
    y = fft(fftshift(vec_))/N
    mag = abs(y)
    
    points = mag[:l]
    ind = arange(l)
    omega = np.sum(points*ind)/np.sum(points)
    start=0
    if omega>1:
        start=1
    delta = mean(angle(y)[start:3])
    return omega,delta


def testEstimator(est,trials=100,noise = False):
    """Test an estimator of omega and delta."""
    t = linspace(-pi,pi,128+1)[:-1]

    oms=[]
    ompreds=[]
    
    dels = []
    delpreds = []

    for i in range(int(trials)):
        omega = 0.5+rand()
        oms.append(omega)
        delta = pi*(rand()-0.5)
        dels.append(delta)
        
        v = cos(omega*t+delta)
        if noise:
            v += 0.1*randn(128)
        om, de = est(v)
        
        ompreds.append(om)
        delpreds.append(de)
        
    oms = array(oms)
    dels = array(dels)
    ompreds = array(ompreds)
    delpreds = array(delpreds)
    omerr = mean(abs(oms-ompreds))
    delerr = mean(abs(dels-delpreds))
    print("MAE for omega: {:.4f}\tMAE for delta: {:.4f}".format(omerr,delerr))
    return omerr,delerr

print("Mean absolute errors for different window sizes:")
print("Estimator 1 without noise:")
for k in 1+arange(10):
    print("k={}".format(k),end='\t')
    testEstimator(lambda v:estimateWD1(v,l=k),trials=1e3, noise=False)

print("Mean absolute errors for different window sizes:")
print("Estimator 1 with noise:")
for k in 1+arange(10):
    print("k={}".format(k),end='\t')
    testEstimator(lambda v:estimateWD1(v,l=k),trials=1e3, noise=True)


def estimateWD2(vec,l=6,window=True):
    """Estimate the value of omega and delta assuming vec contains
    128 samples of cos(omega*t + delta) in (-pi,pi)."""
    
    N = 128
    delta_t = 2*pi/N
    w_max = pi/delta_t
    delta_w = 2*w_max/N
    
    w = linspace(-w_max,w_max,N+1)[:-1]
    
    if window:
        vec_ = vec*fftshift(hamming(arange(N)))
    else:
        vec_ = vec
        
    y = fft(vec_)/N
    mag = abs(y)
    
    points = mag[:l]
    ind = arange(l)
    omega = np.sum(points*ind)/np.sum(points)
    
    
    t = linspace(-pi,pi,N+1)[:-1]
    A = vstack((cos(omega*t),sin(omega*t))).T
    #print(vec.shape)
    a,b = lstsq(A,vec)[0]
    delta = -arctan(b/a)
    
    return omega,delta

print("Estimator 2 without noise:")
for k in 1+arange(10):
    print("k={}".format(k),end='\t')
    testEstimator(lambda v:estimateWD2(v,l=k),trials=1e3, noise=False)

print("Estimator 2 with noise:")
for k in 1+arange(10):
    print("k={}".format(k),end='\t')
    testEstimator(lambda v:estimateWD2(v,l=k),trials=1e3, noise=True)

def chirp(t):
    return cos(16*t*(1.5 + t/(2*pi)))

t = linspace(-pi,pi,1024+1)[:-1]
plot(t,chirp(t))
#scatter(t,chirp(t))
title("Chirp signal")
xlabel("$t$")
ylabel("$y$")
grid()
show()

w,y = plotFFT(chirp, t_range=(-pi,pi), points=1024, wlim=(-80,80), unwrap_=False, plot=True,window=False)

w,y = plotFFT(chirp, t_range=(-pi,pi), points=1024, wlim=(-80,80), unwrap_=False, plot=True,window=True,
             func_name="windowed chirp")

N = 1024
window = 64
n_wins = int(N/window)
delta_t = 2*pi/N
w_max = pi/delta_t
delta_w = 2*w_max/N
    
t = linspace(-pi,pi,N+1)[:-1]
y = chirp(t)

ys=[]
T=26
for i in arange(0,N-window):
    y_ = y[i:i+window]
    Y = 1/window * fftshift(fft(y_))
    ys.append(Y[T:-T])
    
t = linspace(-pi,pi,N)[int(window/2):-int(window/2)]
w = linspace(-w_max,w_max,window+1)[:-1][T:-T]
tt,ww = meshgrid(t,w)
ys = array(ys)


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig2 = figure()
ax2 = fig2.gca(projection='3d')
ax2.set_title("Surface plot of spectrogram of chirp(without hamming window)")
#surf = ax2.plot_surface(tt[:,T:-T],ww[:,T:-T],abs(ys.T)[:,T:-T],cmap=cm.rainbow,
#                       linewidth=0, antialiased=False,rstride=1,cstride=1)
surf = ax2.plot_surface(tt,ww,abs(ys.T),cmap=cm.rainbow,
                       linewidth=0, antialiased=False,rstride=1,cstride=1)

ax2.zaxis.set_major_locator(LinearLocator(10))
ax2.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax2.set_xlabel("$t$")
ax2.set_ylabel("Frequency in rad/s")
ax2.set_zlabel("Magnitude of windowed DFT")
fig2.colorbar(surf, shrink=0.5, aspect=5)
ax2.set_ylim(-80,80)
ax2.view_init(45,15)
    
show()

contourf(tt,ww,abs(ys.T))
grid()
ylim(-80,80)
xlabel("$t$")
ylabel("Frequency in rad/s")
title("Magnitude Spectrogram of chirp(without hamming window)")
show()

contourf(tt,ww,angle(ys.T))
grid()
#ylim(-80,80)
xlabel("$t$")
ylabel("Frequency in rad/s")
title("Phase Spectrogram of chirp(without hamming window)")
show()

N = 1024
window = 64
n_wins = int(N/window)
delta_t = 2*pi/N
w_max = pi/delta_t
delta_w = 2*w_max/N
    
t = linspace(-pi,pi,N+1)[:-1]
y = chirp(t)

ys=[]
T=26
for i in arange(0,N-window):
    y_ = y[i:i+window]*fftshift(hamming(arange(window)))
    Y = 1/window * fftshift(fft(y_))
    ys.append(Y[T:-T])
    
t = linspace(-pi,pi,N)[int(window/2):-int(window/2)]
w = linspace(-w_max,w_max,window+1)[:-1][T:-T]
tt,ww = meshgrid(t,w)
ys = array(ys)


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig2 = figure()
ax2 = fig2.gca(projection='3d')
ax2.set_title("Surface plot of spectrogram of chirp(with hamming window)")
#surf = ax2.plot_surface(tt[:,T:-T],ww[:,T:-T],abs(ys.T)[:,T:-T],cmap=cm.rainbow,
#                       linewidth=0, antialiased=False,rstride=1,cstride=1)
surf = ax2.plot_surface(tt,ww,abs(ys.T),cmap=cm.rainbow,
                       linewidth=0, antialiased=False,rstride=1,cstride=1)

ax2.zaxis.set_major_locator(LinearLocator(10))
ax2.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax2.set_xlabel("$t$")
ax2.set_ylabel("Frequency in rad/s")
ax2.set_zlabel("Magnitude of windowed DFT")
fig2.colorbar(surf, shrink=0.5, aspect=5)
ax2.set_ylim(-80,80)
ax2.view_init(45, 15)
    
show()

contourf(tt,ww,abs(ys.T))
grid()
ylim(-80,80)
xlabel("$t$")
ylabel("Frequency in rad/s")
title("Magnitude Spectrogram of chirp(with hamming window)")
show()

contourf(tt,ww,angle(ys.T))
grid()
#ylim(-80,80)
xlabel("$t$")
ylabel("Frequency in rad/s")
title("Phase Spectrogram of chirp(with hamming window)")
show()

contourf(tt,ww,abs(ys.T))
grid()
#ylim(-80,80)
xlabel("$t$")
ylabel("Frequency in rad/s")
title("Magnitude Spectrogram of chirp(without hamming window)")
show()

contourf(tt,ww,abs(ys.T))
grid()
ylim(-80,80)
xlabel("$t$")
ylabel("Frequency in rad/s")
title("Magnitude Spectrogram of chirp(without hamming window)")
show()

N = 1024
window = 64
n_wins = int(N/window)
delta_t = 2*pi/N
w_max = pi/delta_t
delta_w = 2*w_max/N
    
t = linspace(-pi,pi,N+1)[:-1]
y = chirp(t)

ys=[]
T=0
for i in arange(0,N-window):
    y_ = y[i:i+window]*fftshift(hamming(arange(window)))
    Y = 1/window * fftshift(fft(y_))
    ys.append(Y[T:-T])
    
t = linspace(-pi,pi,N)[int(window/2):-int(window/2)]
w = linspace(-w_max,w_max,window+1)[:-1][T:-T]
tt,ww = meshgrid(t,w)
ys = array(ys)

N = 1024
window = 64
n_wins = int(N/window)
delta_t = 2*pi/N
w_max = pi/delta_t
delta_w = 2*w_max/N
    
t = linspace(-pi,pi,N+1)[:-1]
y = chirp(t)

ys=[]
T=0
for i in arange(0,N-window):
    y_ = y[i:i+window]*fftshift(hamming(arange(window)))
    Y = 1/window * fftshift(fft(y_))
    ys.append(Y)
    
t = linspace(-pi,pi,N)[int(window/2):-int(window/2)]
w = linspace(-w_max,w_max,window+1)[:-1][T:-T]
tt,ww = meshgrid(t,w)
ys = array(ys)

N = 1024
window = 64
n_wins = int(N/window)
delta_t = 2*pi/N
w_max = pi/delta_t
delta_w = 2*w_max/N
    
t = linspace(-pi,pi,N+1)[:-1]
y = chirp(t)

ys=[]
T=0
for i in arange(0,N-window):
    y_ = y[i:i+window]#*fftshift(hamming(arange(window)))
    Y = 1/window * fftshift(fft(y_))
    ys.append(Y)
    
t = linspace(-pi,pi,N)[int(window/2):-int(window/2)]
w = linspace(-w_max,w_max,window+1)[:-1][T:-T]
tt,ww = meshgrid(t,w)
ys = array(ys)

N = 1024
window = 64
n_wins = int(N/window)
delta_t = 2*pi/N
w_max = pi/delta_t
delta_w = 2*w_max/N
    
t = linspace(-pi,pi,N+1)[:-1]
y = chirp(t)

ys=[]
T=0
for i in arange(0,N-window):
    y_ = y[i:i+window]#*fftshift(hamming(arange(window)))
    Y = 1/window * fftshift(fft(y_))
    ys.append(Y)
    
t = linspace(-pi,pi,N)[int(window/2):-int(window/2)]
w = linspace(-w_max,w_max,window+1)[:-1][T:-T]
tt,ww = meshgrid(t,w)
ys = array(ys)

N = 1024
window = 64
n_wins = int(N/window)
delta_t = 2*pi/N
w_max = pi/delta_t
delta_w = 2*w_max/N
    
t = linspace(-pi,pi,N+1)[:-1]
y = chirp(t)

ys=[]
T=0
for i in arange(0,N-window):
    y_ = y[i:i+window]*fftshift(hamming(arange(window)))
    Y = 1/window * fftshift(fft(y_))
    ys.append(Y)
    
t = linspace(-pi,pi,N)[int(window/2):-int(window/2)]
w = linspace(-w_max,w_max,window+1)[:-1][T:-T]
tt,ww = meshgrid(t,w)
ys = array(ys)

w,y = plotFFT(chirp, t_range=(-pi,pi), points=1024, wlim=None, unwrap_=False, plot=True,window=True,
             func_name="windowed chirp")

w,y = plotFFT(chirp, t_range=(-pi,pi), points=1024, wlim=(-80,80), unwrap_=False, plot=True,window=True,
             func_name="windowed chirp")

N = 1024
window = 64
n_wins = int(N/window)
delta_t = 2*pi/N
w_max = pi/delta_t
delta_w = 2*w_max/N
    
t = linspace(-pi,pi,N+1)[:-1]
y = sin(t)

ys=[]
T=0
for i in arange(0,N-window):
    y_ = y[i:i+window]*fftshift(hamming(arange(window)))
    Y = 1/window * fftshift(fft(y_))
    ys.append(Y)
    
t = linspace(-pi,pi,N)[int(window/2):-int(window/2)]
w = linspace(-w_max,w_max,window+1)[:-1][T:-T]
tt,ww = meshgrid(t,w)
ys = array(ys)

N = 1024
window = 64
n_wins = int(N/window)
delta_t = 2*pi/N
w_max = pi/delta_t
delta_w = 2*w_max/N
    
t = linspace(-pi,pi,N+1)[:-1]
y = chirp(t)

ys=[]
T=0
for i in arange(0,N-window):
    y_ = y[i:i+window]*fftshift(hamming(arange(window)))
    Y = 1/window * fftshift(fft(y_))
    ys.append(Y)
    
t = linspace(-pi,pi,N)[int(window/2):-int(window/2)]
w = linspace(-w_max,w_max,window+1)[:-1][T:-T]
tt,ww = meshgrid(t,w)
ys = array(ys)