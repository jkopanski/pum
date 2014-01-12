import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *
from pum.net import *
mu  = 1
eps = 2.56
b   = 2.8 * const.milli 
Z0e = 60
Z0o = 40

modke = Z0e / ( 29.976 * const.pi * np.sqrt( mu / eps))
qe = np.exp( - const.pi * modke)
ke = np.sqrt( qe) * ( ( alg.n_fun( qe) / alg.d_fun( qe)) ** 2)
modko = Z0o / ( 29.976 * const.pi * np.sqrt( mu / eps))
qo = np.exp( - const.pi * modko)
ko = np.sqrt( qo) * ( ( alg.n_fun( qo) / alg.d_fun( qo)) ** 2)

modtest = 2.578092113348 / 1.612441348720
qtest = np.exp( -const.pi * modtest)
ktest = np.sqrt( qtest) * ( ( alg.n_fun( qtest) / alg.d_fun( qtest)) ** 2)
print 'ktest**2 = {}' .format( ktest**2)
print 'ke = {}; ko = {}' .format( ke, ko)

w = ( 2 * b / const.pi) * np.arctanh( np.sqrt( ke * ko))
s = ( 2 * b / const.pi) * np.arctanh( np.sqrt( ke / ko)) - w

print 'w = {} mm; s = {} mm' .format( w / const.milli, s / const.milli)
print '(Z0e, Z0o) = {}' .format( stripline_coupled( w, s, b, 0, mu, eps))
