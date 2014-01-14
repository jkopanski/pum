import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *
from pum.net import *

Z01 = 30.0
Z02 = 75.0
WFS = 1.12
a   = 7 * const.milli
f1  = 2.0 * const.giga
f2  = 3.0 * const.giga
mu  = 1
epsilon = 1
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

print 'V[0] = {0}; V[1] = {1}' .format( np.sqrt( V[0]), V[1])

Z1 = Z01 * np.sqrt( V[0])
Z2 = Z1 * V[1]

print 'Z1 = {0}; Z2 = {1}' .format( Z1, Z2)

b01 = a * np.exp( - Z01 / ( 59.952 *np.sqrt(mu/epsilon)))
b1  = a * np.exp( - Z1 / ( 59.952 *np.sqrt(mu/epsilon)))
b2  = a * np.exp( - Z2 / ( 59.952 *np.sqrt(mu/epsilon)))
b02 = a * np.exp( - Z02 / ( 59.952 *np.sqrt(mu/epsilon)))

print 'b01 = {0}' .format( b01 / const.milli)
print 'b1  = {0}' .format( b1  / const.milli)
print 'b2  = {0}' .format( b2  / const.milli)
print 'b02 = {0}' .format( b02 / const.milli)

print 'l = {0}' .format( const.c / ( np.sqrt( mu * epsilon) * 4 * f0))

R = Z01 / Z02
if R <= 1:
    R = 1 / r
x = f2 / f1

def T21( a2, f):
    theta = const.c / f
    C = ( R - 1) * ( ( R ** 2) + 1) * ( ( R + 1) ** 2) / \
        ( 16 * ( R ** 2) * np.sqrt( R))
    D = 2 * ( R ** 2 + 1) / ( ( R + 1) ** 2)
    E = ( ( R + 1) ** 2)
    F = ( ( R - 1) ** 2) / ( R ** 2 + 1)
    return C * ( \
                 D - ( 2 * np.cos( 2 * a2 * theta)) + \
                 E * np.cos( 2 * theta + ( 2 * a2 * theta)) - \
                 F * ( ( 2 * np.cos( 2 * theta)) - \
                       np.cos( 2 * a2 * theta - ( 2 * theta))))

def equ( p):
    a2, f, L = p
    return ( T21( a2, f1) - L, T21( a2, f) + L, T21( a2, f3) - L)

r = R - 1.5
f1 = 0.629575 * np.exp( -0.115156 * r + \
                        ( 0.004939 * ( r ** 2)) - \
                        ( 0.000074 * ( r ** 3)))
f2 = 0.105558 * np.exp( -0.046644 * r - \
                        ( 0.001213 * ( r ** 2)) + \
                        ( 0.000267 * ( r ** 3)))
f3 = 1.614779 * np.exp( -0.079409 * r + \
                        ( 0.003701 * ( r ** 2)) - \
                        ( 0.000075 * ( r ** 3)))
f4 = 0.251327 - 0.123151 * np.exp( -0.219819 * r + \
                                   ( 0.016291 * ( r ** 2)) - \
                                   ( 0.000646 * ( r ** 3)))
v4 = f1 + ( f2 * ( x - 2))

ts  = v4 / ( 1 + x)
ast = ( f3 + ( f4 * ( 2 - x))) / v4

print 'theta = {0}; a2 = {1}' .format( ts, ast)

a2 = ast
freq = f0
L = T21( ast, f1)
print 'equ = {}' .format( equ(( ast, freq, L)))
print 'a2 = {}; f = {}, L = {}' .format( a2, freq, L)

while np.sqrt( T21( a2, f1)**2 + T21( a2, freq)**2 + T21( a2, f2)**2) > 1e-6:
    df, da, dL = opt.fsolve( equ, ( ast, freq, L))
    a2 += da
    freq += df
    L += dL
    print 'a2 = {}; f = {}, L = {}' .format( da, df, dL)

print 'L1 = {}; L2 = {}; L3 = {}' .format( T21( a2, f1), T21( a2, freq), T21( a2, f2))
