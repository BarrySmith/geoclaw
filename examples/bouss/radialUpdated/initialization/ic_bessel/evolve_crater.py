
from pylab import *
import radial_waves_hankel as RWH

if 1:
    RC = 1500.            # inner radius
    DC = 1000.               # depth of crater
    lip = False

if 0:
    RC = 6676./2             # inner radius
    DC = 2225.               # depth of crater
    lip = False
    
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
     
xlower = 0.
xupper =  40e3
h0 = 4000.

# Initial eta and u

g = 9.81

L = xupper - xlower
x0 = 15e3
mx_H1 = 1000
dx = L/mx_H1
r = linspace(dx/2, xupper-dx/2, mx_H1)
k = linspace(1e-6,0.02,4000)

print('Computing etahat...')
eta0 = eta0_crater(r,RC,DC,lip)
eta0hat = RWH.Htransform(r,eta0,k)

#omega = lambda k: RWH.omega_airy(k,h0)
omega = lambda k: RWH.omega_madsen(k,h0,1/15.)

print('Computing eta, u ...')
t = 0.
eta,u = RWH.eta_u_radial_H(t,r,k,eta0hat,omega,h0,direction='both')

figure(2)
clf()
plot(r,eta0)
plot(r,eta)

def plot_eta(t):
    eta,u = RWH.eta_u_radial_H(t,r,k,eta0hat,omega,h0,direction='both')
    figure(2); clf()
    plot(r,eta)
    
def compare_omega(t):
    omega = lambda k: RWH.omega_swe(k,h0)
    eta_swe,u_swe = \
        RWH.eta_u_radial_H(t,r,k,eta0hat,omega,h0,direction='both')
    omega = lambda k: RWH.omega_airy(k,h0)
    eta_airy,u_airy = \
        RWH.eta_u_radial_H(t,r,k,eta0hat,omega,h0,direction='both')
    omega = lambda k: RWH.omega_madsen(k,h0,1/15.)
    eta_madsen,u_madsen = \
        RWH.eta_u_radial_H(t,r,k,eta0hat,omega,h0,direction='both')
    #omega = lambda k: RWH.omega_madsen(k,h0,0.)
    #eta_madsen0,u_madsen0 = \
    #    RWH.eta_u_radial_H(t,r,k,eta0hat,omega,h0,direction='both')
    omega = lambda k: RWH.omega_sgn(k,h0,1.159)
    eta_sgn,u_sgn = \
        RWH.eta_u_radial_H(t,r,k,eta0hat,omega,h0,direction='both')
    figure(3); clf()
    plot(r, eta_swe, 'k', label='SWE')
    plot(r, eta_airy, 'r', label='Airy')
    #plot(r, eta_madsen0, 'g', label='Madsen B=0')
    plot(r, eta_madsen, 'b', label='Madsen B=1/15')
    plot(r, eta_sgn, 'c', label='SGN, alpha = 1.159')
    legend()
    grid(True)
    title('Time = %.1f seconds' % t)
    savefig('compare_t%s.png' % str(int(t)).zfill(3))
    
