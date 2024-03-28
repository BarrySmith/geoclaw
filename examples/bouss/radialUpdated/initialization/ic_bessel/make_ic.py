from pylab import *
from scipy.special import jv, hankel1

g = 9.81
h0 = 4000.

if 0:
    RC = 9960/2.             # inner radius
    #RD = RC               # no lip
    RD = np.sqrt(2.)*RC  # with lip (PAIR)
    DC = 3502               # depth of crater

if 0:
    RC = 5800/2.             # inner radius
    RD = RC               # no lip
    #RD = np.sqrt(2.)*RC  # with lip (PAIR)
    DC = 10. #2126.              # depth of crater

if 0:
    RC = 11e3             # inner radius
    RD = RC               # no lip
    DC = 10. #5e3               # depth of crater
    
if 0:
    RC = 6676./2             # inner radius (PAIR)
    RD = np.sqrt(2.)*RC  # with lip (PAIR)
    #RD = RC               # no lip
    #DC = 2.*RC/3.         # depth of crater
    DC = 2225.               # depth of crater (PAIR)

if 1:
    RC = 1500.            # inner radius
    RD = RC               # no lip
    DC = 1000.               # depth of crater
    


def F(k):
    """
    F(k,RC,RD) from Ward-Asphaug equation (4)
    """
    F = 4*pi * RD**2 * DC / (RC**2 * k**2) \
        *(jv(2,k*RD) - k*(RD**2 - RC**2) * jv(1,k*RD) / (2*RD))
    return F


def plot_F():
    figure(1,figsize=(9,6))
    clf()
    kmax = 4*8*pi/RC
    k = linspace(0.01*kmax, kmax, 500)
    plot(k, abs(k*F(k)))
    title('F(k,RC,RD) for RC = %g, RD = %g' % (RD,RD))

def omega_WA(k,h):
    """From Ward-Asphaug, agrees with Airy below"""
    return sqrt(g*k*tanh(k*h))

def omega_airy(k,h):
    return k*sqrt(g*h*tanh(k*h)/(k*h))

def c2(kh,B):
    return (1+B*kh**2)/(1 + (B+1/3.)*kh**2)

def omega_madsen(k,h,B):
    return k*sqrt(g*h*c2(k*h,B))

def omega_sgn(k,h,alpha):
    return k*sqrt(g*h*(1 + (alpha-1)* (k*h)**2 /3)  \
                / (1 + alpha*(k*h)**2 /3))

def omega_swe(k,h):
    return k*sqrt(g*h)

def omega(k):
    #omega = omega_WA(k,h0)
    omega = omega_airy(k,h0)
    #omega = omega_sgn(k,h0,alpha=1.)
    #omega = omega_sgn(k,h0,alpha=1.152)
    #omega = omega_madsen(k,h0,B=1/15.)
    #omega = omega_swe(k,h0)
    return omega

def G(k,r,t):
    """
    Integrand of Ward-Asphaug equation (3)
    """
    if type(k) is not ndarray:
        k = array(k)
    if type(r) is not ndarray:
        r = array(r)
    kk,rr = meshgrid(k,r,indexing='ij')
    G = kk/(2*pi) * F(kk) * jv(0,kk*rr) * cos(omega(kk)*t)
    return G

def usurf(r,t,kmax,dk):
    """
    Integrate G*dk using Ward-Asphaug equation (3)
    Returns usurf array with same shape as r
    """
    #dk = kmax/300.
    #dk = kmax/10000.
    k = arange(dk, kmax+dk, dk)
    Gkr = G(k,r,t)
    usurf = dk * sum(Gkr, axis=0)
    return usurf



if 0:
    r = linspace(0,20e3,1000)
    t = 0.
    kmax = 4*8*pi/RC
    dk = kmax/300.
    eta = -usurf(r,t,kmax,dk)
    figure(2,figsize=(9,6))
    clf()
    fill_between(-r/1e3,-h0, eta, color=[.7,.7,1])
    fill_between(r/1e3,-h0, eta, color=[.7,.7,1])
    plot(-r/1e3, eta, 'b')
    plot(r/1e3, eta, 'b')
    grid(True)
    title('Surface eta from usurf(r,t) at t = %.1f' % t)
    xlabel('r (km)')
    ylabel('eta (m)')

def make_fig4():
    kmax = 4*8*pi/RC
    dk = kmax/300.
    r = linspace(0,20e3,1000)
    tt = arange(0,175,10)

    fig,axs = subplots(6,3,figsize=(11,9))
    for jj,t in enumerate(tt):
        i = mod(jj,6)
        j = int((jj-i)/6)
        print('computing eta for ',i,j)
        eta = -usurf(r,t,kmax,dk)
        ax = axs[i,j]
        ax.fill_between(-r/1e3,-h0, eta, color=[.7,.7,1])
        ax.fill_between(r/1e3,-h0, eta, color=[.7,.7,1])
        ax.plot(-r/1e3, eta, 'b')
        ax.plot(r/1e3, eta, 'b')
        ax.text(-18,-3000,'%.0f s\n%.1f m' % (t,abs(eta).max()))
        ax.grid(True)
        ax.set_xlim(-20,20)
        #ax.set_ylim(-h0,h0)
        ax.set_ylim(-20,20)
        if i==5:
            ax.set_xlabel('r (km)')
        else:
            ax.set_xticks([])
        if j==0:
            ax.set_ylabel('eta (m)')
        else:
            ax.set_yticks([])
    


def ploteta_fig5(t):
    r = linspace(0,100e3,2000)
    kmax = 16*8*pi/RC
    dk = kmax / 10000
    eta = -usurf(r,t,kmax,dk)
    figure(2,figsize=(9,6))
    clf()
    plot(-r/1e3, eta, 'b')
    plot(r/1e3, eta, 'b')
    grid(True)

def Fu(k):
    """
    F(k) multiplied by factor giving depth-averaged horizontal velocity
    for each wave number k.
    Does this work for other omega than omega_airy??
    """
    #ufactor = 1/(omega_airy(k,h0)) * (1 - 1/cosh(k*h0))
    ufactor = omega_airy(k,h0) / (k*h0)
    return ufactor*F(k)


def Gu(k,r,t):
    """
    Integrand of Ward-Asphaug equation (3), using Fu in place of F.
    """
    if type(k) is not ndarray:
        k = array(k)
    if type(r) is not ndarray:
        r = array(r)
    kk,rr = meshgrid(k,r,indexing='ij')
    Gu = kk/(2*pi) * Fu(kk) * jv(0,kk*rr) * cos(omega_airy(kk,h0)*t)
    return Gu

def uu(r,t,kmax,dk):
    """
    Integrate Gu*dk using Ward-Asphaug equation (3)
    Returns uu array with same shape as r
    """
    #dk = kmax/300.
    #dk = kmax/10000.
    k = arange(dk, kmax+dk, dk)
    Gukr = Gu(k,r,t)
    uu = dk * sum(Gukr, axis=0)
    return uu


def plot_eta_u(t):
    r = linspace(0,100e3,2000)
    kmax = 8*8*pi/RC
    dk = kmax / 4000
    eta = -usurf(r,t,kmax,dk)
    u = -uu(r,t,kmax,dk)
    u_swe = sqrt(g/h0) * eta
    a = u.max() / u_swe.max()

    figure(2,figsize=(9,6))
    clf()
    xlimits = (-50,50)
    subplot(211)
    plot(-r/1e3, eta, 'b')
    plot(r/1e3, eta, 'b')
    xlim(xlimits)
    grid(True)
    title('Surface eta at t = %.0f sec' % t)

    subplot(212)
    if 1:
        plot(-r/1e3, u, 'b',label='Airy potential')
        plot(r/1e3, u, 'b')

    if 0:
        plot(-r/1e3, u/a, 'b',label='Airy potential / %.1f' % a)
        plot(r/1e3, u/a, 'b')

    plot(-r/1e3, u_swe, 'r',label='SWE evector')
    plot(r/1e3, u_swe, 'r')
    xlim(xlimits)
    legend()
    grid(True)
    title('Depth-averaged u')
    tight_layout()

def plot_Shat():
    kh0 = linspace(-30,30,1000)
    Shat_airy = lambda kh0: sqrt(tanh(kh0)/(kh0))
    Shat_SGN1 = lambda kh0: sqrt(1/(1+kh0**2/3))
    alpha = 1.159
    Shat_SGNa = lambda kh0: sqrt((1+(alpha-1)*kh0**2/3)/(1+alpha*kh0**2/3))
    figure(6,figsize=(11,5))
    clf()
    plot(kh0, Shat_airy(kh0), label='Airy')
    plot(kh0, Shat_SGN1(kh0), label='SGN, alpha=1')
    plot(kh0, Shat_SGNa(kh0), label='SGN, alpha=1.159')
    title('Scale factor for depth-averaged velocity')
    xlabel('k*h0 = depth / wavelength')
    ylabel('Scale factor Shat(k*h0)')
    xlim(-30,30)
    grid(True)
    legend(loc='upper right',framealpha=1)

def Gf(r,fr,k):
    """
    Integrand of Hankel transform of f(r).
    Or inverse transform with (k,fk,r)
    """
    if type(k) is not ndarray:
        k = array(k)
    if type(r) is not ndarray:
        r = array(r)
        fr = array(fr)
    rfr = r*fr
    kk,rr = meshgrid(k,r,indexing='ij')
    kk,rfrf = meshgrid(k,rfr,indexing='ij')
    Gfkr = rfrf * jv(0,kk*rr)
    return Gfkr

def Htransform(r,fr,k):
    """
    Hankel transform if r,fr=f(r) provided, for wave numbers k
    Or inverse transform if (k,fk,r) provided.
    """
    dr = r[1]-r[0]
    Gfkr = Gf(r,fr,k)
    Hf = dr * sum(Gfkr, axis=1)
    return Hf

def Hf(r,fr,k):
    """
    Integrand of unilateral Hankel transform of f(r).
    Or inverse transform with (k,fk,r)
    """
    if type(k) is not ndarray:
        k = array(k)
    if type(r) is not ndarray:
        r = array(r)
        fr = array(fr)
    rfr = r*fr
    kk,rr = meshgrid(k,r,indexing='ij')
    kk,rfrf = meshgrid(k,rfr,indexing='ij')
    Hfkr = rfrf * hankel1(0,kk*rr)
    return Hfkr

def H1transform_old(r,fr,k):
    """
    Unilateral Hankel transform if r,fr=f(r) provided, for wave numbers k
    Or inverse transform if (k,fk,r) provided.
    """
    dr = r[1]-r[0]
    Hfkr = Hf(r,fr,k)
    H1f = dr * sum(Hfkr, axis=1)
    return H1f

def H1transform(r,fr,k,vhankel=0):
    """
    Unilateral Hankel transform if r,fr=f(r) provided, for wave numbers k
    Or inverse transform if (k,fk,r) provided.
    """
    if type(k) is not ndarray:
        # in case only one output value requested
        k = array(k)
    dr = r[1]-r[0]
    rfr = r*fr
    kk,rr = meshgrid(k,r,indexing='ij')
    kk,rfrf = meshgrid(k,rfr,indexing='ij')
    Hfkr = rfrf * hankel1(vhankel,kk*rr)
    H1f = dr * sum(Hfkr, axis=1)
    return H1f
    
def eta_outgoing_radial(t,r,k,etahat1,omega):
    """
    Requires etahat1 = H1transform(r,eta0,k)
        where eta0 is the initial eta evaluated on array r.
    r and k should not include origin
    For in-going, replace -sin term by +sin.
    """
    etahat2 = real(etahat1) * cos(omega(k)*t) - imag(etahat1) * sin(omega(k)*t)
    eta2 = H1transform(k,etahat2,r)
    return eta2

def eta_u_radial_old(t,r,k,etahat,omega,h0,direction='outgoing'):
    etahat2 = zeros(k.shape)
    uhat2 = zeros(k.shape)
    etathat2 = zeros(k.shape)
    if direction in ['outgoing', 'both']:
        etahat2 += real(etahat * exp(1j*omega(k)*t))
        uhat2 += real(etahat * exp(1j*omega(k)*t)) * omega(k) / (h0*k)
        etathat2 += real(1j*omega(k)*etahat * exp(1j*omega(k)*t))
    if direction in ['ingoing', 'both']:
        etahat2 += real(etahat * exp(-1j*omega(k)*t))
        uhat2 += -real(etahat * exp(-1j*omega(k)*t)) * omega(k) / (h0*k)
        etathat2 += real(-1j*omega(k)*etahat * exp(-1j*omega(k)*t))
    if direction == 'both':
        etahat2 *= 0.5
        etathat2 *= 0.5
        uhat2 *= 0.5
    eta2 = H1transform(k,etahat2,r)
    etat2 = H1transform(k,etathat2,r)
    u2 = H1transform(k,uhat2,r)
    
    dr = r[1] - r[0]
    u3 = zeros(u2.shape)
    for j in range(1,len(u3)):
        u3[j] = u3[j-1] + r[j]*etat2[j]
    u3 = -dr/(h0*r)*u3
    #import pdb; pdb.set_trace()
    
    return eta2, u2, etat2, u3
        
def eta_u_radial_H(t,r,k,etahat,omega,h0,direction='outgoing',vhankel=0):
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
    eta = H1transform(k,real(etahat2),r,vhankel)
    etat = H1transform(k,real(etathat2),r,vhankel)
    eta = real(eta)
    etat = real(etat)
    
    dr = r[1] - r[0]
    ru = zeros(eta.shape)
    for j in range(1,len(ru)):
        ru[j] = ru[j-1] - dr*r[j]*etat[j]
    u = ru/(h0*r)
    
    return eta, u
    
def JYtransform(r,fr,k,Bfcn=jv,vhankel=0):
    """
    Hankel transform using Bfcn = jv or yv
    if r,fr=f(r) provided, for wave numbers k
    Or inverse transform if (k,fk,r) provided.
    """
    if type(k) is not ndarray:
        # in case only one output value requested
        k = array(k)
    dr = r[1]-r[0]
    rfr = r*fr
    kk,rr = meshgrid(k,r,indexing='ij')
    kk,rfrf = meshgrid(k,rfr,indexing='ij')
    Jfkr = rfrf * Bfcn(vhankel,kk*rr)
    Jf = dr * sum(Jfkr, axis=1)
    return Jf

def eta2_u2_radial(t,r,k,etahat,omega,h0,vhankel=0):
    
    from scipy.special import jv,yv
    
    omegak = omega(k)
    
    etahat_cos = etahat * cos(omegak*t)
    eta_cosJ = JYtransform(k, etahat_cos, r, jv, vhankel)
    etat_cosY = JYtransform(k, omegak*etahat_cos, r, yv, vhankel)

    etahat_sin = etahat * sin(omegak*t)
    eta_sinY = JYtransform(k, etahat_sin, r, yv, vhankel)
    etat_sinJ = JYtransform(k, -omegak*etahat_sin, r, jv, vhankel)
    
    eta_out = eta_cosJ - eta_sinY
    eta_in  = eta_cosJ + eta_sinY
    etat_out = -etat_sinJ - etat_cosY
    etat_in  = -etat_sinJ + etat_cosY
    
    dr = r[1] - r[0]
    ru_out = zeros(eta_out.shape)
    ru_in = zeros(eta_in.shape)
    for j in range(1,len(ru_out)):
        ru_out[j] = ru_out[j-1] - dr*r[j]*etat_out[j]
        ru_in[j] = ru_in[j-1] - dr*r[j]*etat_in[j]
    u_out = ru_out/(h0*r)
    u_in = ru_in/(h0*r)
    
    return eta_out, u_out, eta_in, u_in
