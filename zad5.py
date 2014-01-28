import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *
from pum.fdm import *

mu  = 1
eps = 1
b   = 8 * const.milli 
t   = 4 * const.milli

Z0 = square_coax( b, t, mu, eps)
print 'Z0(elip) = {}' .format( Z0)
a = b
strip = struct( a, b, 100, 100, mu, eps)
strip.add_plane( 0.0, 0.0, 0.0, 2*b, 0.0)
strip.add_plane( 0.0,   b, 2*a,   b, 0.0)
strip.add_plane(   a,   b,   a,  -b, 0.0)
strip.add_plane(  -a, 0.0,   a, 0.0, 0.0)

strip.add_rect( t, t, a / 2, b / 2, 1.0)
strip.init()
#strip.plot()
print 'k1e-7 = {}' .format( strip.liebmann_quater( 1e-7))
strip.plot()
print 'Z0 = {}' .format( strip.impedance())
