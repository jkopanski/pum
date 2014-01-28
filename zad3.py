import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *

b = 9 * const.milli
eps = 2.04
mu = 1

def find_z( d, mu, eps):
    return cylindrical_flat( d, b, mu, eps) - 30

my_d = alg.newton_raphson( find_z, 0.01 * const.milli, 3.99 * const.milli, tol=1e-10, args = ( 1, 1))
d = opt.newton( find_z, my_d, args = ( 1, 1));

print 'py d = {} mm' .format( d / const.milli)
print 'my d = {} mm' .format( my_d / const.milli)

z0 = cylindrical_flat( my_d, b, 1, 1)
z1 = cylindrical_flat( my_d, b, mu, eps)
print 'Z(mu_r = 1; eps_r = 1) = {}' .format( z0)
print 'Z(mu_r = 2.04; eps_r = 1) = {}' .format( z1)
print 'Z change = {}' .format(z1 - z0)
