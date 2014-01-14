import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *
from pum.net import *
mu  = 1
eps = 2.56
b   = 2.8 * const.milli 
C   = 13
Z0  = 50
f0  = 1.34 * const.giga

k = 10 ** ( - np.abs(C) / 20)
print 'k = {}' .format( k)
Z0e = Z0 * np.sqrt( ( 1 + k) / ( 1 - k))
Z0o = Z0 * np.sqrt( ( 1 - k) / ( 1 + k))
print '(Z0e, Z0o) = {}' .format( Z0e, Z0o)


modke = Z0e / ( 29.976 * const.pi * np.sqrt( mu / eps))
qe = np.exp( - const.pi * modke)
ke = np.sqrt( qe) * ( ( alg.n_fun( qe) / alg.d_fun( qe)) ** 2)
modko = Z0o / ( 29.976 * const.pi * np.sqrt( mu / eps))
qo = np.exp( - const.pi * modko)
ko = np.sqrt( qo) * ( ( alg.n_fun( qo) / alg.d_fun( qo)) ** 2)

w = ( 2 * b / const.pi) * np.arctanh( np.sqrt( ke * ko))
s = ( 2 * b / const.pi) * np.arctanh( np.sqrt( ke / ko)) - w

lamb = const.c / ( np.sqrt(eps) * f0)
print 'lambda = {}; lambda/4 = {}' .format( lamb, lamb / 4)
print 'w = {} mm; s = {} mm' .format( w / const.milli, s / const.milli)
print '(Z0e, Z0o) = {}' .format( stripline_coupled( w, s, b, 0, mu, eps))
