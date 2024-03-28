
from pylab import *
from scipy.special import jv, yv
import make_ic as M

h0 = 4000.

#omega = lambda k: M.omega_airy(k,h0)
#omega = lambda k: M.omega_swe(k,h0)
omega = lambda k: M.omega_madsen(k,h0,1./15.)

#k = linspace(1e-6,0.004,2000)
#r = linspace(1e-3,40e3,1000)
#eta0 = exp(-((r-8e3)/2e3)**2)

k = linspace(1e-6,0.01,4000)
mx_H1 = 1000
#k = linspace(1e-6,0.01,6000)
#mx_H1 = 3000
xupper = 100e3
dx = xupper/mx_H1
r = linspace(dx/2, xupper-dx/2, mx_H1)
#r = linspace(1e-3,100e3,1000)
eta0 = exp(-((r-15e3)/5e3)**2)


figure(1, figsize=(6,7))
clf()

if 0:
    etahat = M.JYtransform(r,eta0,k,jv,0)
    for t in range(0,110,20):
        eta_out,u_out,eta_in,u_in = M.eta2_u2_radial(t,r,k,etahat,omega,h0)
        eta = 0.5*(eta_out + eta_in)
        subplot(211)
        plot(r,eta,label='t = %i' % t)
        subplot(212)
        plot(r,eta_out,label='t = %i' % t)

etahatH = M.H1transform(r,eta0,k)
for t in range(0,210,40):
    eta2,u2 = M.eta_u_radial_H(t,r,k,etahatH,omega,h0,direction='outgoing')
    subplot(211)
    plot(r,eta2,label='H: t = %i' % t)
    subplot(212)
    plot(r,u2,label='H: t = %i' % t)

subplot(211)
title('Surface eta')
legend()
subplot(212)
title('Radial velocity')
legend()
tight_layout()
