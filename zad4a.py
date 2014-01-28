import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *
from pum.fdm import *

mu  = 1
eps = 2.56
b   = 2.8 * const.milli 
Z0  = 50
t   = 0.15 * const.milli

kkp = Z0 / ( 29.976 * const.pi * np.sqrt( mu / eps))
modk = 1 / kkp
q = np.exp( - const.pi * modk)
k = np.sqrt( q) * ( ( alg.n_fun( q) / alg.d_fun( q)) ** 2)
w = ( 2 * b / const.pi) * np.log( ( 1 / k) + np.sqrt( ( 1 / k ** 2) - 1))
a = 10 * b
print 'w = {} mm' .format( w / const.milli)

strip = struct( a, b, 500, 50, mu, eps)
strip.add_plane( 0.0, 0.0, 0.0, 2*b, 0.0)
strip.add_plane( 0.0,   b, 2*a,   b, 0.0)
strip.add_plane(   a,   b,   a,  -b, 0.0)
strip.add_plane(  -a, 0.0,   a, 0.0, 0.0)

strip.add_rect( w, t, a / 2, b / 2, 1.0)
strip.init()
#strip.plot()
print 'k1e-7 = {}' .format( strip.liebmann_quater( 1e-7))
print 'Z0 = {}' .format( strip.impedance())
strip.plot()
