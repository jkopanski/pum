import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *

eps = 2.55
mu  = 1
h   = 1 * const.milli
t   = 0.03 * const.milli
f   = 1 * const.mega
w   = 2.762 * const.milli

print 'Z0 = {}' .format( microstrip( w, t, h, f, mu, eps))

eps = 9.60
mu  = 1
h   = 1 * const.milli
t   = 0.03 * const.milli
f   = 10 * const.giga
w   = 2.201 * const.milli

print 'Z0 = {}' .format( microstrip( w, t, h, f, mu, eps))

def find_w( w):
    return microstrip( w, 0.0035 * const.milli, 1.4 * const.milli, 1.5 * const.giga, 1, 2.56) - 50

my_w = alg.newton_raphson( find_w, 2 * const.milli, 3 * const.milli, (), 1e-10)
#w = opt.newton( find_w, my_w)

#print 'py w = {} mm' .format( w / const.milli)
print 'my w = {} mm' .format( my_w / const.milli)
print 'Z(my w) = {}' .format( microstrip( my_w, 0.0035 * const.milli, 1.4 * const.milli, 1.5 * const.giga, 1, 2.56))

w = my_w
t = 0.0035 * const.milli
h = 1.4 * const.milli
f = 1.5 * const.giga
mu = 1
epsilon = 2.56
u  = w / h
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

print 'lambda = {}' .format( const.c / ( np.sqrt(eps_eff) * f))
