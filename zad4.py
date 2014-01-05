import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *
from pum.net import *
mu  = 1
eps = 2.56
b   = 2.8 * const.milli 

modk = 50 / ( 29.976 * const.pi * np.sqrt( mu / eps))
q = np.exp( - const.pi * modk)
k = np.sqrt( q) * ( ( alg.n_fun( q) / alg.d_fun( q) ** 2))
w = ( 2 * b / const.pi) * np.log( ( 1 / k) + np.sqrt( ( 1 / k ** 2) - 1))

print 'w = {} mm' .format( w / const.milli)

line = net( 12 * const.milli, 6 * const.milli, 7, 4, 6, 3)
line.net[0] = 0.5
line.net[6] = 0.5
for i in range(1, 6):
    line.net[0*7 + i] = 1.0
    line.net[3*7 + i] = 0.0
for j in range(1, 3):
    line.net[j*7 + 0] = 0.0
    line.net[j*7 + 6] = 0.0
print '{}' .format(line.net)
line.init()
print '{}' .format(line.net)

print 'k = {}' .format( line.liebmann( 1e-7))
