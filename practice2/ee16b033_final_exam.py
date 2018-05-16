
from pylab import *

rcParams["figure.figsize"] = 18,9
rcParams["font.size"] = 18
rcParams["text.usetex"] = True

x = linspace(0,10,101)
xx,yy = meshgrid(x,x)
dx=0.1
dy=0.1

sigma = 1e-3
s = zeros((101,101))
s[where(((xx-9)**2+(yy-9)**2<0.8**2)&((xx-9)**2+(yy-9)**2>0.4**2))] = sigma

contourf(xx,yy,s)
ylim(10,0)
show()
Q=300
q = zeros((101,101))
q[where((abs(xx-5)<1)&(abs(yy-5)<1))] = Q
contourf(xx,yy,q)
ylim(10,0)
grid()
show()

def iterate(T,sigarray):
    #print(np.max(sigarray))
    T[1:-1,1:-1] = (0.25*(T[1:-1,0:-2]+
                        T[1:-1,2:]+
                        T[0:-2,1:-1]+
                        T[2:,1:-1])+q[1:-1,1:-1]*dx*dy)/(1+sigarray[1:-1,1:-1]*dx*dy)
    
def bndry(T):
    T[:,0]=T[:,1]
    T[:,-1]=T[:,-2]
    T[0,:]=T[1,:]
    T[-1,:]=T[-2,:]
    #T[where(((xx-9)**2+(yy-9)**2<0.8**2)&((xx-9)**2+(yy-9)**2>0.4**2))]=0
T = zeros((101,101))
oldT = T.copy()
errors = []
N=40000
for i in arange(int(N)):
    iterate(T,s)
    bndry(T)
    T[where(T<0)]=0
    error = np.max(abs(T-oldT))
    errors.append(error)
    oldT = T.copy()

contourf(xx,yy,T,100)
colorbar()
ylim(10,0)
semilogy()
plot(errors)
grid()
NN=10000
x = arange(N)[NN:]
y = log(array(errors))[NN:]
A = c_[x,ones(N)[NN:]]
m,b = lstsq(A,y)[0]
print(m,b)

plot(x,m*x+b)
plot(x,y)
grid()

steps=(log(1e-2)-b)/m
print(steps)

def solve(sig):
    sigarray = zeros((101,101))
    sigarray[where(((xx-9)**2+(yy-9)**2<0.8**2)&((xx-9)**2+(yy-9)**2>0.4**2))] = sig
    print(np.max(sigarray))
    Temp = zeros((101,101))
    oldTemp = Temp.copy()
    errors = []
    N=20000
    for i in arange(int(N)):
        iterate(Temp,sigarray)
        bndry(Temp)
        Temp[where(Temp<0)]=0
        error = np.max(abs(Temp-oldTemp))
        errors.append(error)
        oldT = Temp.copy()
    return Temp,array(errors)
minT = []
maxT = []
for sig in linspace(1e-3,1e5,25):
    Temp,errs = solve(sig)
    minT.append(np.min(Temp))
    maxT.append(np.max(Temp))
    print(sig,minT,maxT)





contourf(xx,yy,T,100)
colorbar()
ylim(10,0)
print(np.max(T))

plot(maxT)

plot(minT)

def solve(sig):
    sigarray = zeros((101,101))
    sigarray[where(((xx-9)**2+(yy-9)**2<0.8**2)&((xx-9)**2+(yy-9)**2>0.4**2))] = sig
    print(np.max(sigarray))
    Temp = zeros((101,101))
    oldTemp = Temp.copy()
    errors = []
    N=40000
    for i in arange(int(N)):
        iterate(Temp,sigarray)
        bndry(Temp)
        Temp[where(Temp<0)]=0
        error = np.max(abs(Temp-oldTemp))
        errors.append(error)
        oldT = Temp.copy()
    return Temp,array(errors)
minT = []
maxT = []
for sig in linspace(1e-3,1e5,25):
    Temp,errs = solve(sig)
    minT.append(np.min(Temp))
    maxT.append(np.max(Temp))
    print(sig,minT,maxT)