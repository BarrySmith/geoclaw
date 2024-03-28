
from pylab import *
from scipy import fft

g = 9.81

Shat_airy = lambda xi,h0: where(abs(xi)>0, sqrt(tanh(xi*h0)/(xi*h0)), 1)
omega_airy = lambda xi,h0: sqrt(g*h0) * xi * Shat_airy(xi,h0)

Shat_sgn = lambda xi,h0,alpha:  sqrt((1 + (alpha-1)* (xi*h0)**2 /3) \
            / (1 + alpha*(xi*h0)**2 /3))
omega_sgn = lambda xi,h0,alpha: sqrt(g*h0) * xi * Shat_sgn(xi,h0,alpha)

Shat_ms = lambda xi,h0,B1: sqrt((1 + B1*(h0*xi)**2) \
            / (1 + (B1 + 1/3.)*(h0*xi)**2))
omega_ms = lambda xi,h0,B1: sqrt(g*h0) * xi * Shat_ms(xi,h0,B1)

# Note: in all casees, omega = sqrt(g*h0) * xi * Shat 
#       and uhat = sqrt(g/h0) * eta_hat * Shat

def eta_u(t, eta0, h0, L, dx, x0, eqn='airy', make_plot=True):
    if eqn=='Airy':
        Shat = lambda xi: Shat_airy(xi,h0)
        omega = lambda xi: sqrt(g*h0) * xi * Shat(xi)
    elif eqn=='SGN1':
        Shat = lambda xi: Shat_sgn(xi,h0,alpha=1)
        omega = lambda xi: sqrt(g*h0) * xi * Shat(xi)
    elif eqn=='SGNa':
        Shat = lambda xi: Shat_sgn(xi,h0,alpha=1.159)
        omega = lambda xi: sqrt(g*h0) * xi * Shat(xi)
    else:
        raise ValueError('Unrecognized eqn')

    x = arange(0, L, dx)
    N = len(x)
    print('N = %i' % N)
    eta = eta0(x) + eta0(x-L)  # make periodic
    etahat = fft.fft(eta)
    xik = fft.fftfreq(N, d=dx) * 2*pi
    shift = exp(-1j*omega(xik)*t)
    etahat_t = etahat * shift
    Shat_k = Shat(xik)
    uhat_t = etahat_t * Shat_k * sqrt(g/h0)

    eta_t = fft.ifft(etahat_t)
    u_t = fft.ifft(uhat_t)

    k0 = argmax(abs(etahat))
    xi0 = xik[k0]
    print('k0 = %i, xi0 = %g' % (k0,xi0))
    cg = (omega(xi0+1e-5)-omega(xi0-1e-5))/2e-5
    print('group velocity for xi0 = %.2f' % cg)

    if make_plot:
    
        figure(13,figsize=(11,8))
        clf()
        subplot(211)
        plot(x-x0,eta_t, 'b', label='eta(t)')
        eta_max = abs(eta).max() * 1.2
        eta_lim = (-eta_max, eta_max)
        plot([cg*t, cg*t],eta_lim, 'r--',label='cg*t')
        ylim(eta_lim)
        legend(loc='lower right',framealpha=1,fontsize=10)
        grid(True)
        ylabel('surface eta (m)')
        title('solution at time t = %.1f using %s' % (t,eqn))

        subplot(212)
        u_swe = eta_t*sqrt(g/h0)
        plot(x-x0, u_swe, 'r', label='SWE eigenvector')
        plot(x-x0,u_t, 'b', label='u(t)')
        legend(loc='lower right',framealpha=1,fontsize=10)
        ylim(sqrt(g/h0)*eta_lim[0], sqrt(g/h0)*eta_lim[1])
        grid(True)
        ylabel('depth-averaged u (m/s)')

    return x, eta_t, u_t


def case1(t):
    """try t = 0, 2000, 4000, 6000"""
    h0 = 4000.
    L = 400e3
    dx = 40.
    x0 = 40e3
    eta0 = lambda x: exp(-((x-x0)/10e3)**2)*cos(2*pi*(x-x0)/5e3)
    x,eta,u = eta_u(t, eta0, h0, L, dx, x0, eqn='Airy')


def case2(t):
    """try t = 0, 500, 1000, 1500"""
    h0 = 4000.
    L = 400e3
    dx = 40.
    x0 = 40e3
    eta0 = lambda x: exp(-((x-x0)/2e3)**2)
    #x,eta,u = eta_u(t, eta0, h0, L, dx, x0, eqn='Airy')
    x,eta,u = eta_u(t, eta0, h0, L, dx, x0, eqn='SGN1')

def case3(t):
    """try t = 0, 500, 1000, 1500"""
    h0 = 4000.
    L = 400e3
    dx = 40.
    x0 = 40e3
    eta0 = lambda x: exp(-((x-x0)/5e3)**2)
    #x,eta,u = eta_u(t, eta0, h0, L, dx, x0, eqn='Airy')
    x,eta,u = eta_u(t, eta0, h0, L, dx, x0, eqn='SGN1')


def case_crater1(t):
    h0 = 4000.
    L = 400e3
    dx = 40.
    x0 = 40e3

    # Ward-Asphaug figure 5
    RC = 5800/2.             # inner radius
    #RD = RC               # no lip
    RD = np.sqrt(2.)*RC  # with lip (PAIR)
    DC = 2126.              # depth of crater
    
    # cut eta0 in half for right-going part:
    eta0 = lambda x: 0.5*where(abs(x-x0)<RD, -DC*(1-(x-x0)**2/RC**2), 0)
    x,eta,u = eta_u(t, eta0, h0, L, dx, x0, eqn='Airy')
 
def compare1(t):
    h0 = 4000.
    L = 400e3
    dx = 40.
    x0 = 40e3
    eta0 = lambda x: exp(-((x-x0)/5e3)**2)
    x,eta1,u1 = eta_u(t, eta0, h0, L, dx, x0, eqn='Airy', make_plot=False)
    x,eta2,u2 = eta_u(t, eta0, h0, L, dx, x0, eqn='SGN1', make_plot=False)
    x,eta3,u3 = eta_u(t, eta0, h0, L, dx, x0, eqn='SGNa', make_plot=False)
    figure(14,figsize=(10,5))
    clf()
    plot(x-x0,eta1,'r',label='Airy')
    plot(x-x0,eta2,'b',label='SGN alpha=1.')
    plot(x-x0,eta3,'k',label='SGN alpha=1.159')
    legend()
    grid(True)
    
    
