import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *

a = 7    * const.milli
b = 3.04 * const.milli
c = 0

print 'Z coax: {}' .format( coax_z( a, b, 1, 1))
print 'Z skew: {}' .format( skew_coax_z( c, a, b, 1, 1))

print 'Z error {}' .format( (coax_z( a, b, 1, 1) - skew_coax_z( c, a, b, 1, 1)) / coax_z( a, b, 1, 1) * 100)

def analytical( a, b, mu_r, epsilon_r):
    d = coax_z(a, b, mu_r, epsilon_r) - 5
    k = 59.952 * np.sqrt(mu_r/epsilon_r)
    return np.sqrt( (( a*a) - ( 2*a*b*np.cosh(d/k)) + ( b*b)) / 4)

def find_coax( c):
    return skew_coax_z( c, a, b, 1, 1) - coax_z( a, b, 1, 1) + 5

c = opt.newton( find_coax, 0.8 * const.milli);
my_c = alg.newton_raphson( find_coax, 0, 2 * const.milli, tol=1e-10)

print 'analytical shift {} mm' .format( analytical( a, b, 1, 1) / const.milli)
print 'numerical shift {} mm' .format( c / const.milli)
print 'Z skewed: {}' .format( skew_coax_z( c, a, b, 1, 1))
print 'shift my newt: {} mm' .format( my_c / const.milli)
print 'my error: {}' .format( np.abs( my_c - c))
