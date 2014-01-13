import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *
from pum.net import *

Z01 = 25.0
Z02 = 75.0
WFS = 1.12
f1  = 2.0 * const.giga
f2  = 3.0 * const.giga

f0    = ( f2 + f1) / 2
w     = ( ( f2 - f1)) / f0
gamma = ( WFS - 1) / ( WFS + 1)
r     = Z01 / Z02
if r <= 1:
    r = 1 / r
print 'r = {0}, gamma = {1}' .format( r, gamma)

nx = np.arccosh( ( r - 1) / ( gamma * ( r + 1))) / \
     np.arccosh( 1 / np.cos( const.pi * ( 1 - ( w / 2)) / 2))
n  = np.ceil( nx)

print 'np = {0}; n = {1}' .format( nx, n)

if n == 1:
    V = np.sqrt( r)
elif n == 2:
    mu   = np.sin( const.pi * w / 4)
    c    = ( ( r - 1) * ( mu ** 2)) / ( 2 * ( 2 - ( mu ** 2)))
    V1 = np.sqrt( ( c ** 2) + r) + c
    V2 = r / V1
    V  = (V1, V2)

print 'V[0] = {0}; V[1] = {1}' .format( V[0], V[1])

Z1 = Z01 * np.sqrt( V[0])
Z2 = Z1 * V[1]

print 'Z1 = {0}; Z2 = {1}' .format( Z1, Z2)
