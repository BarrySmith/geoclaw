"""
Utilities for solving radial depth-averaged linear Airy, shallow water,
and Boussinesq-type equations for a given dispersion relation omega(k).
Implementation of Hankel transforms using direct integration.
"""


from pylab import *
from scipy.special import hankel1

g = 9.81

# Define several dispersion relations. To use one of these in eta_u_radial_H,
# define a function omega(k) using, e.g.
#   h0 = 4000.  # constant ocean depth
#   omega = lambda k: omega_airy(k,h0)


def omega_airy(k,h):
    """
    Airy solution dispersion relation,
    omega(k) = sqrt(g*k*tanh(k*h)) = k*sqrt(g*h*tanh(k*h)/(k*h))
    """
    return sqrt(g*k*tanh(k*h))
    
def omega_madsen(k,h,B):
    """
    Madsen-Sorenson dispersion relation with parameter B.
    For B=0 this agrees with SGN for alpha=1.
    """
    c2 = (1+B*(k*h)**2)/(1 + (B+1/3.)*(k*h)**2)
    return k*sqrt(g*h*c2)

def omega_sgn(k,h,alpha):
    """
    Serre Green Naghdi dispersion relation with parameter alpha.
    Classic case has alpha=1.
    """
    return k*sqrt(g*h*(1 + (alpha-1)* (k*h)**2 /3)  \
                / (1 + alpha*(k*h)**2 /3))

def omega_swe(k,h):
    """
    Linearized shallow water equation dispersion relation.
    """
    return k*sqrt(g*h)
    
    
def Htransform(r,fr,k,vhankel=0):
    """
    Hankel transform if r and fr=f(r) provided, for wave numbers k.
    The same formulas provide the inverse transform if (k,fk,r) are provided
    instead of (r,fr,k).
    """
    if type(k) is not ndarray:
        # in case only one output value requested
        k = array(k)
    dr = r[1]-r[0]
    rfr = r*fr
    kk,rr = meshgrid(k,r,indexing='ij')
    kk,rfrf = meshgrid(k,rfr,indexing='ij')
    Hfkr = rfrf * hankel1(vhankel,kk*rr)
    Hf = dr * sum(Hfkr, axis=1)
    return Hf
    

def eta_u_radial_H(t,r,k,etahat,omega,h0,direction='outgoing'):
    """
    Compute eta,u for the radially symmetric linear equation determined by
    the dispersion relation omega(k) provided as a function of one variable,
    on an ocean with constant depth h0.
    
    t is the time at which eta, and u are desired.

    etahat is the Hankel transform of the initial data eta(r,0), and is
    an array of the same shape as k.
        
    If direction=='both', the initial data for u is u(r,0)=0.
    If direction=='outgoing' or 'ingoing', then u(r,0) is chosen so that the 
    wave is "purely" outgoing or ingoing, with minimal energy going in the
    other direction.
    
    Note that the Hankel transform of eta_t(r,t) is also computed and
    transformed back to eta_t(r,t), and then u(r,t) is computed
    by integrating the mass equation in radial coordinates,
        eta_t + (1/r)*(ru)_r = 0.
    """
    etahat2 = zeros(k.shape,dtype=complex128)
    etathat2 = zeros(k.shape,dtype=complex128)
    if direction in ['outgoing', 'both']:
        etahat2 += etahat * exp(1j*omega(k)*t)
        etathat2 += 1j*omega(k)*etahat * exp(1j*omega(k)*t)
    if direction in ['ingoing', 'both']:
        etahat2 += etahat * exp(-1j*omega(k)*t)
        etathat2 += -1j*omega(k)*etahat * exp(-1j*omega(k)*t)
    if direction == 'both':
        etahat2 *= 0.5
        etathat2 *= 0.5
    eta = Htransform(k,real(etahat2),r,vhankel=0)
    etat = Htransform(k,real(etathat2),r,vhankel=0)
    eta = real(eta)
    etat = real(etat)
    
    dr = r[1] - r[0]
    ru = zeros(eta.shape)
    for j in range(1,len(ru)):
        ru[j] = ru[j-1] - dr*r[j]*etat[j]
    u = ru/(h0*r)
    
    return eta, u
