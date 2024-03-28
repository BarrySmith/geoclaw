
from pylab import *
from scipy import fft

L = 100e3
dx = 40.
x = arange(dx/2, L, dx)
N = len(x)
print('N = %i' % N)

xi0 = 0.002
beta = 10e3

h0 = 1000.
Shat = lambda xi: where(xi>0, sqrt(tanh(xi*h0)/(xi*h0)), 1)
print('x = %.3f, Shat(xi0) = %.3f' % (xi0,Shat(xi0)))

u = cos(xi0*x) * exp(-(x/beta)**2)
#u = cos(xi0*x) * exp(-((x-25e3)/beta)**2)
u = cos((x/4e3)**2) * exp(-(x/30e3)**2)
uhat = fft.dct(u)
xik = arange(0,N,1.)*pi/L
Shat_k = Shat(xik)
Suhat = Shat_k * uhat
Su = fft.idct(Suhat)
ratio = Su.max()/u.max()
print('Ratio Su_max / u_max = %.3f' % ratio)

figure(11,figsize=(7,7))
clf()
subplot(211)
plot(x,u, 'b', label='original u')
plot(x,Su, 'r', label='after scaling\nfor h0 = %g' % h0)
legend(loc='upper right',framealpha=1)
wavelength = 2*pi/xi0
kmax = xi0*L/pi
grid(True)
#title('xi0 = %.3f, wavelength = %.1f, kmax = %.1f' % (xi0,wavelength,kmax))

subplot(212)
plot(uhat, 'b', label='u hat')
plot(Suhat, 'r', label='after scaling\nfor h0 = %g' % h0)
legend(loc='upper right',framealpha=1)
#xlim(0,3*kmax)
grid(True)


