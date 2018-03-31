
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
    
from sympy import *
import scipy.signal as sp

H,s=symbols('H(s) s')
init_printing()

def lowpass(R1=10e3,R2=10e3,C1=1e-9,C2=1e-9,G=1.586,Vi=1):
    """Solve the given lowpass filter circuit for a given input Vi."""
    A=Matrix([[0,0,1,-1/G],
              [-1/(1+s*R2*C2),1,0,0],
              [0,-G,G,1],
              [-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    
    b=Matrix([0,0,0,-Vi/R1])
    
    V=A.inv()*b
    return (A,b,V)

def bodePlot(H_s,w_range=(0,8),points=800):
    """Plot the magnitude and phase of H_s over the given range of frequencies."""
    
    w = logspace(*w_range,points)
    h_s = lambdify(s,H_s,'numpy')
    H_jw = h_s(1j*w)
    
    # find mag and phase
    mag = 20*np.log10(np.abs(H_jw))
    phase = angle(H_jw,deg = True)
    
    eqn = Eq(H,simplify(H_s))
    display(eqn)
    
    fig,axes = plt.subplots(1,2,figsize=(18,6))
    ax1,ax2 = axes[0],axes[1]
    
    # mag plot
    ax1.set_xscale('log')
    ax1.set_ylabel('Magntiude in dB')
    ax1.set_xlabel('$\omega$ in rad/s')
    ax1.plot(w,mag)
    ax1.grid()
    ax1.set_title("Magnitude of $H(j \omega)$")
    
    # phase plot
    ax2.set_ylabel('Phase in degrees')
    ax2.set_xlabel('$\omega$ in rad/s')
    ax2.set_xscale('log')
    ax2.plot(w,phase)
    ax2.grid()
    ax2.set_title("Phase of $H(j \omega)$")
    
    plt.show()

A,b,V=lowpass()
Vo = V[3]
bodePlot(Vo)

def symToTransferFn(Y_s):
    """
    Convert a sympy rational polynomial into one that can be used for scipy.signal.lti.
    
    Returns a tuple (num,den) which contains the coefficients of s of
    the numerator and denominator polynomials.
    """
    
    # get the polynomial coefficients after simplification
    Y_ss = expand(simplify(Y_s))
    n,d = fraction(Y_ss)
    n,d = Poly(n,s), Poly(d,s)
    num,den = n.all_coeffs(), d.all_coeffs()
    num,den = [float(f) for f in num], [float(f) for f in den]

    return num,den

def findQ(H_s):
    """Find the quality factor of the input transfer function assuming it is second order."""
    nl,dl = symToTransferFn(H_s)
    syst = sp.lti(nl,dl)
    p1,p2 = syst.poles[0], syst.poles[1]
    return np.sqrt(abs(p1*p2))/abs(p1+p2)

print("Q factor of low pass filter: {:.4f}".format(findQ(lowpass()[-1][-1])))

#G1 = symbols('G1')
def highpass(R1=10e3,R3=10e3,C1=1e-9,C2=1e-9,G=1.586,Vi=1):
    """Solve the given highpass filter circuit for a given input Vi."""
    
    A=Matrix([[0,-1,0,1/G],
        [s*C2*R3/(s*C2*R3+1),0,-1,0],
        [0,G,-G,1],
        [-s*C2-1/R3-s*C1,1/R3,0,1/R1]])

    b=Matrix([0,0,0,-Vi*s*C1])

    V=A.inv()*b
    #V_lim = [limit(v,G1,oo) for v in V]
    return (A,b,V)

A,b,V = highpass()
V_highpass =simplify(V[3])
bodePlot(V_highpass)

print("Q factor of high pass filter: {:.4f}".format(findQ(highpass()[-1][-1])))

def inverseLaplace(Y_s,t=None):
    """
    Finds the inverse laplace transform of a sympy expression using sp.impulse. 
    """
    
    # load the step response as a system in scipy.signal
    num,den = symToTransferFn(Y_s)
    
    # evaluate in time domain
    t,y = sp.impulse((num,den),T=t)
    return t,y

def plotFilterOutputs(laplace_in=None, time_domain_fn=None, 
                      lp_range=(0,1e-3), hp_range=(0,1e-3), points=1e3,
                     input_name="Input"):
    """
    Plot the time domain outputs of the two active filters to a given input in the 
    laplace domain or the time domain.
    """
    
    t_lp = linspace(*lp_range,points)
    t_hp = linspace(*hp_range,points)

    
    if laplace_in != None:
        
        A,b,V_lowpass = lowpass(Vi=laplace_in)
        t_lp,y_lp = inverseLaplace(V_lowpass[-1],t=t_lp)
        
        
        A,b,V_highpass = highpass(Vi=laplace_in)
        t_hp,y_hp = inverseLaplace(V_highpass[-1],t=t_hp)
        
    elif time_domain_fn != None:
        
        A,b,V_lowpass = lowpass()
        lowsys = symToTransferFn(V_lowpass[-1])
        t_lp,y_lp,svec = sp.lsim(lowsys, time_domain_fn(t_lp), t_lp)
        
         
        A,b,V_highpass = highpass()
        highsys = symToTransferFn(V_highpass[-1])
        t_hp,y_hp,svec = sp.lsim(highsys, time_domain_fn(t_hp), t_hp)
    
    else:
        print("No input given.")
        
     
    fig,axes = plt.subplots(1,2,figsize=(18,6))
    ax1,ax2 = axes[0],axes[1]
    
    # low pass response plot
    ax1.set_ylabel('$V_o$')
    ax1.set_xlabel('$t$')
    ax1.plot(t_lp,y_lp)
    ax1.grid()
    ax1.set_title("Response of low pass filter to {}".format(input_name))
    
    # high pass response plot
    ax2.set_ylabel('$V_o$')
    ax2.set_xlabel('$t$')
    ax2.plot(t_hp,y_hp)
    ax2.grid()
    ax2.set_title("Response of high pass filter to {}".format(input_name))
   
    plt.show()
    return t_lp,y_lp, t_hp,y_hp

t_lp,y_lp, t_hp,y_hp = plotFilterOutputs(laplace_in=1/s,input_name="unit step")
print("Steady state value of low pass filter step response: {:.4f}".format(y_lp[-1]))
print("Steady state value of high pass filter step response: {:.4f}".format(y_hp[-1]))

def vi1(t):
    """Sum of low frequency and high frequency sinusoids"""
    u_t = 1*(t>0)
    return (np.sin(2000*np.pi*t)+np.cos(2e6*np.pi*t)) * u_t

t = linspace(0,1e-3,1e5)
plt.title("Sum of low and high frequency sinusoids")
plt.xlabel("$t$")
plt.ylabel("$v_i(t)$")
plt.plot(t,vi1(t))
plt.grid()
plt.show()

a = plotFilterOutputs(time_domain_fn=vi1,input_name="sum of sinusoids",points=1e5,hp_range=(0,1e-4),lp_range=(0,1e-5))

a = plotFilterOutputs(time_domain_fn=vi1,input_name="sum of sinusoids",points=1e5,hp_range=(0,1e-5),lp_range=(0,3e-3))

def input_f(t,decay=0.5,freq=1.5):
    """Exponentially decaying cosine function."""
    u_t = 1*(t>0)
    return np.cos(freq*t)*np.exp(-decay*t) * u_t

t = linspace(0,1e-3,1e5)
plt.title("High frequency damped sinusoid")
plt.xlabel("$t$")
plt.ylabel("$v_i(t)$")
plt.plot(t,input_f(t,decay=3e3,freq=1e7))
plt.grid()
plt.show()

a = plotFilterOutputs(time_domain_fn=lambda t:input_f(t,decay=3e3,freq=1e7),points=1e5
                  ,input_name='\nhigh frequency damped sinusoid')

t = linspace(0,1.8e-1,1e5)
plt.title("Low frequency damped sinusoid")
plt.xlabel("$t$")
plt.ylabel("$v_i(t)$")
plt.plot(t,input_f(t,decay=1e1,freq=1e3))
plt.grid()
plt.show()

a=plotFilterOutputs(time_domain_fn=lambda t:input_f(t,decay=1e1,freq=1e3),lp_range=(0,1.8e-1),points=1e5
                  ,input_name='\nlow frequency damped sinusoid')