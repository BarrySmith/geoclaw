
from pylab import *
from scipy import fft

L = 320e3
dx = 20.
x = arange(0, L, dx)
#x = arange(-L/2, L/2, dx)
N = len(x)
print('N = %i' % N)

xi0 = 0.0005
beta = 20e3

xi0 = 0.00001
beta = 5e3

h0 = 4000.
g = 9.81
Shat = lambda xi: where(abs(xi)>0, sqrt(tanh(xi*h0)/(xi*h0)), 1)
omega = lambda xi: sqrt(g*h0) * xi * Shat(xi)
print('xi0 = %.3f, Shat(xi0) = %.3f' % (xi0,Shat(xi0)))

#xt = x
#u = cos(xi0*(xt-0)) * exp(-((xt-0)/beta)**2) \
#   + cos(xi0*(xt-L)) * exp(-((xt-L)/beta)**2)
#u = cos(xi0*x) * exp(-((x-25e3)/beta)**2)
#u = cos((x/4e3)**2) * exp(-(x/30e3)**2)

def u0(x):
    u0a = lambda x: cos(xi0*x) * exp(-(x/beta)**2)
    return u0a(x) + u0a(x-L)  # make periodic

u = u0(x)
    
uhat = fft.fft(u)
#xik = arange(0,N,1.)*2*pi/L
xik = fft.fftfreq(N, d=dx) * 2*pi
Shat_k = Shat(xik)
#Shat_k = 1
Suhat = Shat_k * uhat
Su = fft.ifft(Suhat)
ratio = Su.max()/u.max()
ratio_hat = Suhat.max()/uhat.max()
print('Ratio Suhat_max / uhat_max = %.3f' % ratio)
print('Ratio Su_max / u_max = %.3f' % ratio)

figure(11,figsize=(7,7))
clf()
subplot(211)
plot(x,u, 'b', label='original u')
plot(x,Su, 'r', label='after scaling\nfor h0 = %g' % h0)
legend(loc='upper right',framealpha=1)
wavelength = 2*pi/xi0
kmax = xi0*L/(2*pi)
grid(True)
title('xi0 = %.3f, wavelength = %.1f, kmax = %.1f' % (xi0,wavelength,kmax))

subplot(212)
plot(uhat, 'b', label='u hat')
plot(Suhat, 'r', label='after scaling\nfor h0 = %g' % h0)
legend(loc='upper right',framealpha=1)
xlim(0,3*kmax)
grid(True)

def plot_u(t):
    
    u = u0(x)
    uhat = fft.fft(u)
    xik = fft.fftfreq(N, d=dx) * 2*pi
    shift = exp(-1j*omega(xik)*t)
    Shat_k = Shat(xik)
    Suhat = Shat_k * uhat * shift
    Su = fft.ifft(Suhat)

    cg = (omega(xi0+1e-5)-omega(xi0-1e-5))/2e-5
    print('group velocity = %.2f' % cg)

    figure(12,figsize=(7,7))
    clf()
    subplot(211)
    #plot(x,u, 'b', label='u(0)')
    plot(x,Su, 'r', label='u(t)')
    plot([cg*t, cg*t],[-1,1], 'k--',label='cg*t')
    legend(loc='upper right',framealpha=1)
    wavelength = 2*pi/xi0
    kmax = xi0*L/(2*pi)
    grid(True)
    title('u at time t = %.1f' % t)

    subplot(212)
    #plot(uhat, 'b', label='u hat')
    plot(Suhat, 'r', label='at time t')
    legend(loc='upper right',framealpha=1)
    xlim(0,3*kmax)
    grid(True)

def plot_eta_u(t, eta0):
    
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

    figure(13,figsize=(11,8))
    clf()
    subplot(211)
    plot(x,eta_t, 'b', label='eta(t)')
    plot([cg*t, cg*t],[-1,1], 'k--',label='cg*t')
    legend(loc='lower right',framealpha=1,fontsize=10)
    grid(True)
    ylabel('surface eta (m)')
    title('solution at time t = %.1f' % t)

    subplot(212)
    u_swe = eta_t*sqrt(g/h0)
    plot(x, u_swe, 'r', label='SWE eigenvector')
    plot(x,u_t, 'b', label='u(t)')
    legend(loc='lower right',framealpha=1,fontsize=10)
    ylim(-sqrt(g/h0), sqrt(g/h0))
    grid(True)
    ylabel('depth-averaged u (m/s)')
