import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *

eps = 1
mu  = 1
h   = 8 * const.milli
x = 0.5483
y = 0.7882
d = x * h
s = y * h

print '(Z0e, Z0o) = {}' .format( cylindrical_flat_coupled( s, d, h, mu, eps))

Z0e = 60.0
Z0o = 40.0

def V( x, y):
    ind = x * h
    ins = y * h
    Z = cylindrical_flat_coupled( ins, ind, h, mu, eps)
    return ( Z[0] - Z0e, Z[1] - Z0o)

Z0 = np.sqrt( Z0e * Z0o)
k = ( Z0e - Z0o) / ( Z0e + Z0o)
x = ( 4 / const.pi) * \
    np.exp( ( -Z0) / \
            ( 59.952 * \
              np.sqrt( ( 0.987) - \
                       ( 0.171 * k) - \
                       ( 1.723 * ( k**3)))))
r = ( 4 / ( const.pi * x)) ** ( 0.001 + ( 1.117 * k) - ( 0.683 * ( k ** 2)))
y = ( 1 / const.pi) * \
    np.log( ( r + 1) / ( r - 1)) - x
d = x * h
s = y * h

print '(Z0e, Z0o) = {}' .format( cylindrical_flat_coupled( s, d, h, mu, eps))

print 'Z0 = {}; k = {}; r = {}; x = {}; y = {}' .format( Z0, k, r, x, y)
delta = 1e-8
(V1, V2) = V( x, y)
while ( V1**2 + V2**2) > ( Z0e * Z0o * 1e-6):
    (V1dxp, V2dxp) = V( x + delta, y)
    (V1dxn, V2dxn) = V( x - delta, y)
    V1dx = ( V1dxp - V1dxn) / ( 2 * delta)
    V2dx = ( V2dxp - V2dxn) / ( 2 * delta)
    (V1dyp, V2dyp) = V( x, y + delta)
    (V1dyn, V2dyn) = V( x, y - delta)
    V1dy = ( V1dyp - V1dyn) / ( 2 * delta)
    V2dy = ( V2dyp - V2dyn) / ( 2 * delta)
    J = V1dx * V2dy - ( V1dy * V2dx)
    if J == 0:
        x += delta
        y -= delta
        continue
    (V1, V2) = V( x, y)
    x -= ( ( V1 * V2dy) - ( V2 * V1dy)) / J
    y += ( ( V1 * V2dx) - ( V2 * V1dx)) / J

d = x * h
s = y * h
print 's = {}; d = {}' .format( s / const.milli, d / const.milli)
print '(Z0e, Z0o) = {}' .format( cylindrical_flat_coupled( s, d, h, mu, eps))
