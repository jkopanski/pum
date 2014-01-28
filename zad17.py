import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
import matplotlib.pyplot as plt
from pum.lines import *
from pum.net import *

Z0  = 50.0
f0  = 2.8 * const.giga
fa  = 3.2 * const.giga
w   = 0.1
Lr  = 0.2
La  = 30.0
eps = 1
mu  = 1

Lap = 10 ** ( La / 10)
Lrp = 10 ** ( Lr / 10)
wa  = 2 * ( fa - f0) / f0

nx = np.arccosh( np.sqrt( ( Lap - 1) / (Lrp - 1))) / \
     np.arccosh( np.sin( const.pi * wa / 4) /\
                 np.sin( const.pi * w  / 4))
n = np.ceil( nx)
print 'nx = {}; n = {}' .format( nx, n)

x = np.log( 1 / np.tanh( Lr / 17.37))
y = np.sinh( x / ( 2 * n))
a = [float for x in range(int(n) + 2)]
b = [float for x in range(int(n) + 2)]
g = [float for x in range(int(n) + 2)]
el = [float for x in range(int(n) + 2)]

for k in range(int(n) + 2):
    a[k] = np.sin( ( 2 * k - 1) * const.pi / ( 2 * n))
    b[k] = ( y ** 2) + ( np.sin( k * const.pi / n) ** 2)
    if k == 0:
        g[k] = 1.0
    elif k == 1:
        g[k] = 2 * a[k] / y
    elif k == n + 1:
        if ( k % 2) == 0:
            g[k] = ( 1 / np.tanh( x / 4)) ** 2
        else:
            g[k] = 1
    else:
        g[k] = 4 * a[ k - 1] * a[ k] / ( b[ k - 1] * g[ k - 1])

for k in range(int(n) + 2):
    print 'g{} = {}' .format( k, g[k])

Y0    = 1 / Z0
theta = ( 1 - ( w / 2)) * const.pi / 2
d     = 0.5
# w1 = 2 * const.pi * f1
y  = [float for x in range(1, int(n) + 2)]
N  = [float for x in range(1, int(n) + 2)]
If = [float for x in range(1, int(n) + 2)]
Ib = [float for x in range(1, int(n) + 2)]
w1p = 1
Ca  = 2 * g[1] * d
If[1] = g[0] * np.sqrt( Ca / g[2])
Ib[1] = If[1]
for i in range( 1, int(n) + 1):
    if n == 1:
        If[1], Ib[1] = g[0] * np.sqrt( Ca / g[i+1])
    else:    
        If[i] = g[0] * Ca / np.sqrt( g[i] * g[i+1])
        Ib[i] = g[0] * np.sqrt( Ca * g[i+1] / ( g[0] * g[i-1]))
    N[i]  = np.sqrt( ( If[i] ** 2) + \
                     ( ( g[0] * w1p * Ca * np.tan( theta) / 2) ** 2))
    if i == n:
        y[i] = Y0 * w1p * ( g[i] * g[i+1] - ( d * g[0] * g[1])) * \
               np.tan( theta) + ( Y0 * ( N[i-1] - Ib[i]))
    elif i != 1:
        y[i]  = Y0 * ( N[i-1] + N[i] - Ib[i] - If[i])
    else:
        y[i] = Y0 * ( g[0] * ( 1 - d) * g[1] * w1p * np.tan( theta) + \
                      N[i] - If[i])
for k in range(1, int(n) + 1):
    print 'Z{} = {}' .format( k, 1/y[k])

def find_w( Z0, b, t, mu, eps):
    m  = 6 * ( b - t) / ( 3 * b - t)
    B  = np.exp( Z0 * np.sqrt( eps) / 30)
    W  = 8 * ( ( b - t) / const.pi) * ( np.sqrt( B + 0.568) / ( B - 1))
    dw = t * \
         ( 1 - \
           ( np.log( ( ( t / ( 2 * b - t)) ** 2) + \
                     ( ( ( 0.0796 * t) / ( W - (0.26*t))) ** m) / 2))) / const.pi
    return W - dw

mu  = 1
eps = 2.56
b   = 2.8 * const.milli 
t   = 0.15 * const.milli

for k in range(1, int(n) + 1):
    print 'w{} = {}' .format( k, find_w( 1/y[k], b, t, mu, eps) / const.milli)

print 'l = {}' .format( const.c / ( np.sqrt(eps) *f0) / 4 / const.milli)

#     if ( k % 2) == 1:
#         el[k] = Z0 * g[k] / w1
#     else:
#         el[k] = g[k] / ( Z0 * w1)

# for k in range( int(n) + 2):
#     print( 'g{} = {}; el{} = {}' .format( k, g[k], k, el[k]))

# d0 = D * np.exp( - Z0 / ( 59.952 *np.sqrt(1/1)))
# dl = D * np.exp( - Zl / ( 59.952 *np.sqrt(mu/eps)))
# dh = D * np.exp( - Zh / ( 59.952 *np.sqrt(1/1)))
# print 'd0 = {}; dl = {}; dh = {}' .format( d0 / const.milli, dl / const.milli, dh / const.milli)

# Cf  = const.pi * D * ( 13 + ( 0.2 * D / dh)) * \
#       ( ( ( dl - dh) / ( D - dh)) ** 2.65) * const.pico
# Cf0 = const.pi * D * ( 13 + ( 0.2 * D / dh)) * \
#       ( ( ( d0 - dh) / ( D - dh)) ** 2.65) * const.pico
# print 'Cf0 = {}; Cf = {}' .format( Cf0, Cf)

# vl = const.c / np.sqrt(eps)
# vh = const.c

# # Pierwsze przyblizenie
# for k in range( 1, int(n) + 1):
#     if ( k % 2) == 1:
#         l[k] = np.arcsin( w1 * el[k] / Zh) * vh / w1
#     else:
#         l[k] = w1 * el[k] * Zl * vl / w1

# for i in range( 10):
#     for k in range( 1, int(n) + 1):
#         if k == 1:
#             l[k] = np.arcsin( ( w1 * el[k] - \
#                                 ( ( Zl * l[ k + 1] * w1) / ( 2 * vl))) / Zh) * vh / w1
#         elif k == ( int(n)):
#             l[k] = np.arcsin( ( w1 * el[k] - \
#                                 ( ( Zl * l[ k - 1] * w1) / ( 2 * vl))) / Zh) * vh / w1
#         else:
#             if ( k % 2) == 1:
#                 l[k] = np.arcsin( ( w1 * el[k] - \
#                                     ( ( Zl * l[ k - 1] * w1) / ( 2 * vl)) - \
#                                     ( ( Zl * l[ k + 1] * w1) / ( 2 * vl))) / Zh) * vh / w1
#             else:
#                 l[k] = ( w1 * el[k] - \
#                          ( 2 * Cf * w1) - \
#                          ( ( l[ k - 1] * w1) / ( 2 * vh * Zh)) - \
#                          ( ( l[ k - 1] * w1) / ( 2 * vh * Zh))) * Zl * vl / w1


# lo = ( ( Z0 / Zh) ** 2) * ( Zh * Cf0 * vh + ( l[1] / 2))
# print 'lo = {}' .format( lo / const.milli)
# for k in range( 1, int(n) + 1):
#     if k == 1:
#         l[k] += lo
#     elif k == 7:
#         l[k] += lo
#     print 'l{} = {}' .format( k, l[k] / const.milli)
