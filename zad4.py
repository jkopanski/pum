import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *
from pum.net import *
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

print 'w = {} mm' .format( w / const.milli)

# line = net( 12 * const.milli, 6 * const.milli, 7, 4, 6, 3)
# line.net[0].u = 0.5
# line.net[6].u = 0.5
# for i in range(1, 6):
#     line.net[0*7 + i].u = 1.0
#     line.net[3*7 + i].u = 0.0
# for j in range(1, 3):
#     line.net[j*7 + 0].u = 0.0
#     line.net[j*7 + 6].u = 0.0
# line.init()

# print 'k = {}' .format( line.liebmann( 1e-7))

#line.plot()

x = 10 * w / 2
y = b / 2
xp = int( np.ceil( x / ( 0.01 * const.milli)))
yp = int( np.ceil( y / ( 0.01 * const.milli)))
wp = int( np.ceil( w / 2 / x * xp))
tp = int( np.ceil( t / 2 / y * yp))
print( 'xp ={}; yp = {}; wp = {}; tp = {}' .format( xp, yp, wp, tp))
zad4 = net( x, y, xp+1, yp+1, wp, tp)
for j in range( yp):
    for i in range( xp):
        zad4.net[0][i] = met_point(0.0)
        zad4.net[j][0] = met_point(0.0)
        if i > ( xp - wp):
            if j > ( yp - tp):
                zad4.net[j][i] = met_point(1.0)

for i in range( xp+1):
    zad4.net[yp][i] = imag_point(zad4.net[yp - 2][i])

for j in range( yp):
    zad4.net[j][xp] = imag_point(zad4.net[j][xp - 2])

zad4.plot()
zad4.init()
zad4.plot()
print 'k = {}' .format( zad4.liebmann( 1e-2))
zad4.plot()
