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
C   = 3.01
Z0  = 50
f0  = 1.34 * const.giga
R   = 1

y1 = np.sqrt( 10 ** ( - C / 10))
y2 = np.sqrt( 1 - ( y1 ** 2))
Z1 = Z0 / y1
Z2 = Z0 / y2

print( 'y1 = {}; y2 = {}; Z1 = {}; Z2 = {}' .format( y1, y2, Z1, Z2))

def find_w1( w):
    return microstrip( w, t, h, f0, mu, eps) - Z1
def find_w2( w):
    return microstrip( w, t, h, f0, mu, eps) - Z2

w1 = alg.newton_raphson( find_w1, 0.5 * const.milli, 2 * const.milli, (), tol=1e-10)
w2 = alg.newton_raphson( find_w2, 0.5 * const.milli, 2 * const.milli, (), 1e-10)

print( 'w1 = {}; w2 = {}' .format( w1 / const.milli, w2 / const.milli))
print( 'spr')
print( 'Z1 = {}; Z2 = {}' \
       .format( microstrip( w1, t, h, f0, mu, eps), \
                microstrip( w2, t, h, f0, mu, eps)))

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
k  = 0.363 * g * np.exp( -4.6 * u)
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
k  = 0.363 * g * np.exp( -4.6 * u)
m  = 0.525 / ( ( 1 + ( 0.157 * fn))**20)
n  = 0.27488 + ( 0.6315 * u) + ( m * u) - ( 0.065683 * np.exp( -8.7513 * u))
p  = n * c * ( ( ( 1.844 * fn) + ( k * d * fn))**1.5763)
eps_0 = ( ( epsilon + 1) / 2) + ( ( ( epsilon - 1) / 2) * ( ( 1 + ( 10 / u))**( -a * b)))
eps_eff = ( eps_0 + ( epsilon * p)) / ( 1 + p)

print 'lambda2 = {}' .format( const.c / ( np.sqrt(eps_eff) * f))
print 'lambda2/4 = {}' .format( const.c / ( np.sqrt(eps_eff) * f) / 4)

def coupling(f, f0, Z0, y1, y2):
    theta  = ( const.pi / 2) * f / f0
    A  = ( y2 ** 2) - \
         ( ( y1 ** 2) * np.tan( 3 * theta / 2) * np.tan( theta / 2))
    Ap = ( y2 ** 2) - \
         ( ( y1 ** 2) * ( 1 / np.tan( 3 * theta / 2)) * \
           ( 1 / np.tan( theta / 2)))
    B  = y1 * y2 * ( 1 / np.tan( theta)) * \
         ( np.tan( 3 * theta / 2) + np.tan( theta / 2))
    Bp = y1 * y2 * ( 1 / np.tan( theta)) * \
         ( ( 1 / np.tan( 3 * theta / 2)) + ( 1 / np.tan( theta / 2)))
    C  = y1 * ( np.tan( 3 * theta / 2) + np.tan( theta / 2))
    Cp = y1 * ( ( 1 / np.tan( 3 * theta / 2)) + ( 1 / np.tan( theta / 2)))
    D  = y1 * ( np.tan( 3 * theta / 2) - np.tan( theta / 2))
    Dp = y1 * ( ( 1 / np.tan( 3 * theta / 2)) - ( 1 / np.tan( theta / 2)))
    E  = -2 * y2 * ( 1 / np.tan( theta))
    Spp = ( 1 - A - B + ( 1j * D)) / \
          ( 1 + A + B + ( 1j * ( C + E)))
    Spm = ( 1 - Ap + Bp - ( 1j * Dp)) / \
          ( 1 + Ap - Bp - ( 1j * ( Cp - E)))
    S34 = ( Spp - Spm) / 2
    return 20 * np.log10( 1 / np.absolute( S34))

freq = np.linspace( 1.25 * const.giga, 1.45 * const.giga)
line, = plt.plot( freq, coupling( freq, f0, Z0, y1, y2))
plt.show()
