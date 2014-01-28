import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
import matplotlib.pyplot as plt
from pum.lines import *
from pum.net import *

mu  = 1
eps = 4.34
h   = 1.4 * const.milli 
t   = 0.035 * const.milli
C   = 3.9
Z0  = 50.0
f0  = 1.34 * const.giga
R   = 1.0

k = 1/ np.sqrt( 10 ** ( C / 10) - 1)
print( 'k = {}' .format( k))
Z1 = Z0 / k
Z2 = Z0 * np.sqrt( R / ( 1 + ( k ** 2)))
Z3 = Z0 * R / k

print( 'Z1 = {}; Z2 = {}; Z3 = {}' .format( Z1, Z2, Z3))

def find_w1( w):
    return microstrip( w, t, h, f0, mu, eps) - Z1
def find_w2( w):
    return microstrip( w, t, h, f0, mu, eps) - Z2
def find_w3( w):
    return microstrip( w, t, h, f0, mu, eps) - Z3

w1 = alg.newton_raphson( find_w1, 1 * const.milli, 3 * const.milli, (), tol=1e-10)
w2 = alg.newton_raphson( find_w2, 3 * const.milli, 5 * const.milli, (), 1e-10)
w3 = alg.newton_raphson( find_w3, 1 * const.milli, 3 * const.milli, (), 1e-10)
print( 'w1 = {}; w2 = {}; w3 = {}' .format( w1 / const.milli, w2 / const.milli, w3 / const.milli))

print( 'spr')
print( 'Z1 = {}; Z2 = {}; Z3 = {}' \
    .format( microstrip( w1, t, h, f0, mu, eps), \
             microstrip( w2, t, h, f0, mu, eps), \
             microstrip( w3, t, h, f0, mu, eps)))

f = f0
epsilon = eps
u  = w1 / h
if t != 0:
    du = ( t / ( 2 * const.pi * h)) * np.log( 1 + ( ( 4 * np.exp(1) * h) / ( t * ( ( 1 / np.tanh( np.sqrt( 6.517 * u)))**2)))) * ( 1 + ( 1 / np.cosh( np.sqrt( epsilon - 1))))
    u += du
a  = 1 + ( ( 1 / 49) * np.log( ( ( u**4) + ( ( u / 52)**2)) / ( ( u**4) + 0.432))) + ( ( 1 / 18.7) * np.log( 1 + ( ( u / 18.1)**3)))
b  = 0.564 * ( ( ( epsilon - 0.9) / ( epsilon + 3))**0.053)
c  = 0.33622 * ( 1 - np.exp( -0.03442 * epsilon))
d  = 3.751 - ( 2.75 * np.exp( - ( ( epsilon / 15.916)**8)))
fn = f * h * ( 10**-7)
g  = 1 - np.exp( - ( ( fn / 3.87)**4.97))
#k  = 0.363 * g * np.exp( -4.6 * u)
m  = 0.525 / ( ( 1 + ( 0.157 * fn))**20)
n  = 0.27488 + ( 0.6315 * u) + ( m * u) - ( 0.065683 * np.exp( -8.7513 * u))
p  = n * c * ( ( ( 1.844 * fn) + ( k * d * fn))**1.5763)
eps_0 = ( ( epsilon + 1) / 2) + ( ( ( epsilon - 1) / 2) * ( ( 1 + ( 10 / u))**( -a * b)))
eps_eff = ( eps_0 + ( epsilon * p)) / ( 1 + p)

print 'lambda1 = {}' .format( const.c / ( np.sqrt(eps_eff) * f))
print 'lambda1/4 = {}' .format( const.c / ( np.sqrt(eps_eff) * f) / 4)

u  = w2 / h
if t != 0:
    du = ( t / ( 2 * const.pi * h)) * np.log( 1 + ( ( 4 * np.exp(1) * h) / ( t * ( ( 1 / np.tanh( np.sqrt( 6.517 * u)))**2)))) * ( 1 + ( 1 / np.cosh( np.sqrt( epsilon - 1))))
    u += du
a  = 1 + ( ( 1 / 49) * np.log( ( ( u**4) + ( ( u / 52)**2)) / ( ( u**4) + 0.432))) + ( ( 1 / 18.7) * np.log( 1 + ( ( u / 18.1)**3)))
b  = 0.564 * ( ( ( epsilon - 0.9) / ( epsilon + 3))**0.053)
c  = 0.33622 * ( 1 - np.exp( -0.03442 * epsilon))
d  = 3.751 - ( 2.75 * np.exp( - ( ( epsilon / 15.916)**8)))
fn = f * h * ( 10**-7)
g  = 1 - np.exp( - ( ( fn / 3.87)**4.97))
#k  = 0.363 * g * np.exp( -4.6 * u)
m  = 0.525 / ( ( 1 + ( 0.157 * fn))**20)
n  = 0.27488 + ( 0.6315 * u) + ( m * u) - ( 0.065683 * np.exp( -8.7513 * u))
p  = n * c * ( ( ( 1.844 * fn) + ( k * d * fn))**1.5763)
eps_0 = ( ( epsilon + 1) / 2) + ( ( ( epsilon - 1) / 2) * ( ( 1 + ( 10 / u))**( -a * b)))
eps_eff = ( eps_0 + ( epsilon * p)) / ( 1 + p)

print 'lambda2 = {}' .format( const.c / ( np.sqrt(eps_eff) * f))
print 'lambda2/4 = {}' .format( const.c / ( np.sqrt(eps_eff) * f) / 4)

u  = w3 / h
if t != 0:
    du = ( t / ( 2 * const.pi * h)) * np.log( 1 + ( ( 4 * np.exp(1) * h) / ( t * ( ( 1 / np.tanh( np.sqrt( 6.517 * u)))**2)))) * ( 1 + ( 1 / np.cosh( np.sqrt( epsilon - 1))))
    u += du
a  = 1 + ( ( 1 / 49) * np.log( ( ( u**4) + ( ( u / 52)**2)) / ( ( u**4) + 0.432))) + ( ( 1 / 18.7) * np.log( 1 + ( ( u / 18.1)**3)))
b  = 0.564 * ( ( ( epsilon - 0.9) / ( epsilon + 3))**0.053)
c  = 0.33622 * ( 1 - np.exp( -0.03442 * epsilon))
d  = 3.751 - ( 2.75 * np.exp( - ( ( epsilon / 15.916)**8)))
fn = f * h * ( 10**-7)
g  = 1 - np.exp( - ( ( fn / 3.87)**4.97))
#k  = 0.363 * g * np.exp( -4.6 * u)
m  = 0.525 / ( ( 1 + ( 0.157 * fn))**20)
n  = 0.27488 + ( 0.6315 * u) + ( m * u) - ( 0.065683 * np.exp( -8.7513 * u))
p  = n * c * ( ( ( 1.844 * fn) + ( k * d * fn))**1.5763)
eps_0 = ( ( epsilon + 1) / 2) + ( ( ( epsilon - 1) / 2) * ( ( 1 + ( 10 / u))**( -a * b)))
eps_eff = ( eps_0 + ( epsilon * p)) / ( 1 + p)

print 'lambda3 = {}' .format( const.c / ( np.sqrt(eps_eff) * f))
print 'lambda3/4 = {}' .format( const.c / ( np.sqrt(eps_eff) * f) / 4)

def coupling(f, f0, Z0, R, k):
    Z1 = Z0 / k
    Z2 = Z0 * np.sqrt( R / ( 1 + ( k ** 2)))
    Z3 = Z0 * R / k
    z2 = Z2 / Z0

    theta  = ( const.pi / 2) * f / f0
    c      = np.cos( theta)
    s      = np.sin( theta)
    b1e    = np.tan( theta / 2) * Z0 / Z1
    b3e    = np.tan( theta / 2) * Z0 / Z3
    b1o    = - Z0 / ( Z1 * np.tan( theta / 2))
    b3o    = - Z0 / ( Z3 * np.tan( theta / 2))

    RDe    = ( c * ( R + 1)) - ( s * z2 * ( b1e + ( b3e * R)))
    RDo    = ( c * ( R + 1)) - ( s * z2 * ( b1o + ( b3o * R)))
    XDe    = ( ( c * R * ( b1e + b3e)) + \
               ( s * z2) + \
               ( ( s * R) / z2) - \
               ( s * b1e * b3e * z2 * R))
    XDo    = ( ( c * R * ( b1o + b3o)) + \
               ( s * z2) + \
               ( ( s * R) / z2) - \
               ( s * b1o * b3o * z2 * R))

    RDe    = RDe / ( 2 * np.sqrt( R))
    RDo    = RDo / ( 2 * np.sqrt( R))
    XDe    = XDe / ( 2 * np.sqrt( R))
    XDo    = XDo / ( 2 * np.sqrt( R))
    modS14 = ( 1.0 / 2.0) * \
             ( np.sqrt( ( ( ( RDo - RDe) ** 2) + ( ( XDo - XDe) ** 2))/ \
                        ( ( (RDe ** 2) + (XDe ** 2)) * \
                        ( (RDo ** 2) + (XDo ** 2)))))
    return 20*np.log10( 1.0 / modS14)

freq = np.linspace( 1 * const.giga, 1.68 * const.giga)
line, = plt.plot( freq, coupling( freq, f0, Z0, R, k))
plt.xlabel( 'czestotliwosc [Hz]')
plt.ylabel( 'sprzezenie C [dB]')
plt.show()
