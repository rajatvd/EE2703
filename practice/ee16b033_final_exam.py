
from pylab import *

rcParams["figure.figsize"] = 18,9
rcParams["font.size"] = 18
rcParams["text.usetex"] = True
N = 100
a = 10
phi_ = linspace(0,2*pi,N+1)[:-1]
r_ = c_[a*cos(phi_),a*sin(phi_),zeros(N)]

I =  c_[a*cos(phi_)*sin(phi_),-a*cos(phi_)*cos(phi_),zeros(N)]

quiver(r_[:,0],r_[:,1],I[:,0],I[:,1],scale=700,color='red')
grid()
xlim(-18,18)
ylim(-12,12)
show()

x = arange(3)-1
z = arange(1,1001,1)
xx,yy,zz = meshgrid(x,x,z)
print(xx.shape)
r = zeros((3,3,1000,3))
r[:,:,:,0]=xx
r[:,:,:,1]=yy
r[:,:,:,2]=zz
print(r.shape)
print(r[2,0,23])


R=norm(tile(r,100).reshape(3,3,1000,100,3)-r_,axis=-1)
R.shape

def calc(l):
    return norm(r-r_[l],axis=-1)

A_x = sum(cos(phi_)*exp(0.1j*R)*I[:,0]/R,axis=-1)
A_y = sum(cos(phi_)*exp(0.1j*R)*I[:,1]/R,axis=-1)

#Bz = (A[2,1,:,1]-A[1,2,:,0]+A[0,1,:,1]-A[1,0,:,0])/4
Bz = (A_y[1,2,:]-A_y[1,0,:]-(A_x[2,1,:]-A_x[0,1,:]))/4
loglog()
grid()
plot(z,abs(Bz))
show()

K=50
y = log(abs(Bz))[K:]
x = log(z[K:])
A = c_[x,ones(1000)[K:]]
m,b = lstsq(A,y)[0]
plot(x,y)
print(m,b)
plot(x,m*x+b)
grid()
show()

import scipy.signal as sp

