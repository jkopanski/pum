import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
import matplotlib.pyplot as plt
from pum.lines import *
from pum.net import *

Z0  = 50.0
Zl  = 10.0
Zh  = 120.0
D   = 7 * const.milli
f1  = 1    * const.giga
fa  = 1.43 * const.giga
Lr  = 0.2
La  = 30.0
eps = 2.05
mu  = 1

Lap = 10 ** ( La / 10)
Lrp = 10 ** ( Lr / 10)

nx = np.arccosh( np.sqrt( ( Lap - 1) / (Lrp - 1))) / np.arccosh( fa / f1)
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
            g[k] = ( 1.0 / np.tanh( x / 4.0)) ** 2
        else:
            g[k] = 1
    else:
        g[k] = 4 * a[ k - 1] * a[ k] / ( b[ k - 1] * g[ k - 1])

w1 = 2 * const.pi * f1

for k in range(1, int(n) + 1):
    if ( k % 2) == 1:
        el[k] = Z0 * g[k] / w1
    else:
        el[k] = g[k] / ( Z0 * w1)

for k in range( int(n) + 2):
    print( 'g{} = {}; el{} = {}' .format( k, g[k], k, el[k]))

d0 = D * np.exp( - Z0 / ( 59.952 *np.sqrt(1/1)))
dl = D * np.exp( - Zl / ( 59.952 *np.sqrt(mu/eps)))
dh = D * np.exp( - Zh / ( 59.952 *np.sqrt(1/1)))
print 'd0 = {}; dl = {}; dh = {}' .format( d0 / const.milli, dl / const.milli, dh / const.milli)

Cf  = const.pi * D * ( 13 + ( 0.2 * D / dh)) * \
      ( ( ( dl - dh) / ( D - dh)) ** 2.65) * const.pico
Cf0 = const.pi * D * ( 13 + ( 0.2 * D / dh)) * \
      ( ( ( d0 - dh) / ( D - dh)) ** 2.65) * const.pico
print 'Cf0 = {}; Cf = {}' .format( Cf0, Cf)

vl = const.c / np.sqrt(eps)
vh = const.c

# Pierwsze przyblizenie
l = [float for x in range(1, int(n) + 2)]
for k in range( 1, int(n) + 1):
    if ( k % 2) == 1:
        l[k] = np.arcsin( w1 * el[k] / Zh) * vh / w1
    else:
        l[k] = w1 * el[k] * Zl * vl / w1

for k in range( 1, int(n) + 1):
    if ( k % 2) == 0:
        l[k] = ( w1 * el[k] - \
                 ( 2 * Cf * w1) - \
                 ( ( l[ k - 1] * w1) / ( 2 * vh * Zh)) - \
                 ( ( l[ k + 1] * w1) / ( 2 * vh * Zh))) * Zl * vl / w1

for i in range( 50):
    for k in range( 1, int(n) + 1):
        if k == 1:
            l[k] = np.arcsin( ( w1 * el[k] - \
                                ( ( Zl * l[ k + 1] * w1) / ( 2 * vl))) / Zh) * vh / w1
        elif k == ( int(n)):
            l[k] = np.arcsin( ( w1 * el[k] - \
                                ( ( Zl * l[ k - 1] * w1) / ( 2 * vl))) / Zh) * vh / w1
        else:
            if ( k % 2) == 1:
                l[k] = np.arcsin( ( w1 * el[k] - \
                                    ( ( Zl * l[ k - 1] * w1) / ( 2 * vl)) - \
                                    ( ( Zl * l[ k + 1] * w1) / ( 2 * vl))) / Zh) * vh / w1
    for k in range( 1, int(n) + 1):
            if ( k % 2) == 0:
                l[k] = ( w1 * el[k] - \
                         ( 2 * Cf * w1) - \
                         ( ( l[ k - 1] * w1) / ( 2 * vh * Zh)) - \
                         ( ( l[ k + 1] * w1) / ( 2 * vh * Zh))) * Zl * vl / w1


lo = ( ( Z0 / Zh) ** 2) * ( Zh * Cf0 * vh + ( l[1] / 2))
print 'lo = {}' .format( lo / const.milli)
for k in range( 1, int(n) + 1):
    print 'l{} = {}' .format( k, l[k] / const.milli)
