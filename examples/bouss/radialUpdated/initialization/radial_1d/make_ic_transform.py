
from pylab import *
from scipy.interpolate import interp1d

import sys
sys.path.insert(0,'../ic_bessel')
from radial_waves_hankel import eta_u_radial_H, Htransform, \
                                omega_madsen, omega_swe, omega_airy

xedges = loadtxt('celledges.txt',skiprows=1)[:,0]
xcell = 0.5*(xedges[:-1] + xedges[1:])

xlower = xedges[0]
xupper =  xedges[-1]
L = xupper - xlower
g = 9.81
h0 = 4000.

# Initial eta and u

omega = lambda k: omega_madsen(k,h0,1./15.)
#omega = lambda k: omega_swe(k,h0)
#omega = lambda k: omega_airy(k,h0)



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

if 1:
    # Gaussian ring:
    x0 = 15e3  # ring
    #x0 = 0.   # Gaussian
    #direction = 'both'
    direction = 'outgoing'
    eta0_ring = lambda x: 1*exp(-((x-x0)/5e3)**2)
    #r = linspace(1e-3,80e3,4000)
    mx_H1 = 1000
    dx = L/mx_H1
    r_H1 = linspace(dx/2, xupper-dx/2, mx_H1)
    r = linspace(dx/2, xupper-dx/2, mx_H1)
    k = linspace(1e-6,0.01,4000)
    eta0 = eta0_ring(r)



print('Computing eta0hat...')
eta0hat = Htransform(r,eta0,k)

def make_data_and_plots():

    figure(10); clf(); plot(k,eta0hat)
    fname = 'etahat.png'
    savefig(fname)
    print('Created ',fname)

    t0 = 0.

    eta,u = eta_u_radial_H(t0,r,k,eta0hat,omega,h0,direction=direction)
    etaimax = abs(imag(eta)).max()
    print('etaimax = %g' % etaimax)

    eta_fcn = interp1d(r,eta,fill_value="extrapolate")
    eta_cell = eta_fcn(xcell)
    u_fcn = interp1d(r,u,fill_value="extrapolate")
    u_cell = u_fcn(xcell)
    etau = vstack((real(eta_cell),real(u_cell))).T

    fname = 'eta_u.txt'
    savetxt(fname,etau,fmt='%.8e  ')
    print('Created ',fname)

    figure(11); clf();
    plot(r,eta0,label='original eta0')
    plot(r,eta,label='eta at t0=%.1f' % t0)
    plot(xcell,eta_cell,label='eta in cells')
    plot(xcell,u_cell*sqrt(h0/g), label='u * sqrt(h0/g)')
    grid(True)
    legend()
    fname = 'eta_u.png'
    savefig(fname)
    print('Created ',fname)

if __name__=='__main__':

    make_data_and_plots()

