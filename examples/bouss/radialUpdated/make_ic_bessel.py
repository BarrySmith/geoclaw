
from pylab import *

import sys,os

try:
    BoussDev = os.environ['BoussDev']
except:
    print("*** Need to define environment variable BoussDev to path")

sys.path.insert(0,BoussDev+'/initialization/ic_bessel')

from radial_waves_hankel import eta_u_radial_H, Htransform, \
                                omega_madsen, omega_swe, omega_sgn

xlower = 0.
xupper =  100e3
h0 = 4000.

# Initial eta and u

g = 9.81

#omega = lambda k: omega_madsen(k,h0,1./15.); print('Using Madsen')
#omega = lambda k: omega_sgn(k,h0,1.); print('Using SGN, alpha=1')
omega = lambda k: omega_sgn(k,h0,1.153); print('Using SGN, alpha=1.153')
#omega = lambda k: omega_swe(k,h0); print('Using SWE')

L = xupper - xlower

def eta0_crater(r,RC,DC,lip):
    """
    Initial conditions for a parabolic crater with radius RC, depth DC
    r is an array of radial distance values.  RC,DC,r are all in meters.
    If lip==True then a lip extends out to RD=sqrt(2)*RC.
    """
    if lip:
        RD = sqrt(2.)*RC
    else:
        RD = RC
    eta = -DC*(1 - (r/RC)**2)
    eta = where(r>RD, 0, eta)
    return eta

if 0:
    # crater
    RC = 1500.            # inner radius
    DC = 1000.            # depth of crater
    lip = False
    direction = 'both'
    mx_H1 = 4000
    dx = L/mx_H1
    r = linspace(dx/2, xupper-dx/2, mx_H1)
    k = linspace(1e-6,0.02,4000)

    eta0 = eta0_crater(r,RC,DC,lip)


if 0:
    x0 = 15e3  # ring
    #x0 = 0.   # Gaussian
    #direction = 'both'
    direction = 'outgoing'
    eta0_ring = lambda x: 1*exp(-((x-x0)/5e3)**2)
    mx_H1 = 1000
    dx = L/mx_H1
    r = linspace(dx/2, xupper-dx/2, mx_H1)
    k = linspace(1e-6,0.01,4000)
    eta0 = eta0_ring(r)

if 1:
    x0 = 0.   # Gaussian
    direction = 'both'
    eta0_gaussian = lambda x: 1*exp(-((x-x0)/2e3)**2)
    mx_H1 = 1000
    dx = L/mx_H1
    r = linspace(dx/2, xupper-dx/2, mx_H1)
    k = linspace(1e-6,0.01,4000)
    eta0 = eta0_gaussian(r)

# ----------------------------

etahat1 = Htransform(r,eta0,k)

def make_data_and_plots():
    figure(10); clf(); 
    plot(k,etahat1)
    title('Hankel transform of eta0')
    grid(True)
    fname = 'etahat.png'
    savefig(fname)
    print('Created ',fname)

    ketahat = vstack((k,etahat1)).T
    fname = 'etahat.data'
    savetxt(fname,ketahat,fmt='%.8e  ')
    print('Created ',fname)


    t = 0.
    eta,u = eta_u_radial_H(t,r,k,etahat1,omega,h0,direction=direction)
    etaimax = abs(imag(eta)).max()
    print('etaimax = %g' % etaimax)
    etau = vstack((r,eta0,real(u))).T
    #etau = vstack((r,real(eta),real(u))).T
    fname = 'starting.data'
    savetxt(fname,etau,fmt='%.8e  ',comments='',header='%i\n%i' % (t,mx_H1))
    print('Created ',fname)

    figure(11); clf(); 
    plot(r,eta0,label='eta0')
    plot(r,eta,label='eta after transforms')
    plot(r,u*sqrt(h0/g), label='u * sqrt(h0/g)')
    grid(True)
    legend()
    xlim(-1e3,20e3)
    fname = 'eta_u.png'
    savefig(fname)
    print('Created ',fname)

if __name__ == '__main__':

    make_data_and_plots()


