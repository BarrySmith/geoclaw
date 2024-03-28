
from pylab import *
from clawpack.geoclaw_1d import nonuniform_grid_tools

import sys
#sys.path.insert(0,'../ic_fft')
#from planewave_fft import eta_u
#sys.path.insert(0,'../ic_bessel')
#from make_ic import eta_u_radial_H, H1transform
#from radial_waves_hankel import eta_u_radial_H, Htransform, \
#                                omega_madsen, omega_swe


xlower = 0.
xupper =  75e3
g = 9.81
h0 = 4000.
xzpairs = [(xlower,-h0),   # left edge
           (xupper,-h0)]  # right edge

topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

mx = 3000
hmin = 10.

nonuniform_grid_tools.make_celledges_cfl(xlower, xupper, mx, topo_fcn,
        hmin, fname='celledges.txt', plot_topo=True)

if 0:
    # Initial eta and u

    omega = lambda k: omega_madsen(k,h0,1./15.)

    L = xupper - xlower
    dx = L/mx


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

        mx_H1 = 1000
        dx = L/mx_H1
        r = linspace(dx/2, xupper-dx/2, mx_H1)
        k = linspace(1e-6,0.02,4000)

        eta0 = eta0_crater(r,RC,DC,lip)

    if 0:
        # Gaussian ring:
        x0 = 15e3  # ring
        #x0 = 0.   # Gaussian
        #direction = 'both'
        direction = 'outgoing'
        eta0_ring = lambda x: 1000*exp(-((x-x0)/5e3)**2)
        #r = linspace(1e-3,80e3,4000)
        mx_H1 = 2000
        dx = L/mx_H1
        r_H1 = linspace(dx/2, xupper-dx/2, mx_H1)
        r = linspace(dx/2, xupper-dx/2, mx)
        k = linspace(1e-6,0.01,4000)
        eta0 = eta0_ring(r)



    print('Computing eta0hat...')
    eta0hat = Htransform(r,eta0,k)


    figure(10); clf(); plot(k,eta0hat)
    fname = 'etahat.png'
    savefig(fname)
    print('Created ',fname)

    t = 0.
    eta,u = eta_u_radial_H(t,r,k,eta0hat,omega,h0,direction=direction)
    etaimax = abs(imag(eta)).max()
    print('etaimax = %g' % etaimax)
    etau = vstack((real(eta),real(u))).T
    fname = 'eta_u.txt'
    savetxt(fname,etau,fmt='%.8e  ')
    print('Created ',fname)

    figure(11); clf();
    plot(r,eta,label='eta')
    plot(r,u*sqrt(h0/g), label='u * sqrt(h0/g)')
    grid(True)
    legend()
    fname = 'eta_u.png'
    savefig(fname)
    print('Created ',fname)
