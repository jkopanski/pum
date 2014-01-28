import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
import matplotlib.pyplot as plt
from pum.lines import *
from pum.net import *

Z0  = 50.0
b   = 8 * const.milli
t   = 0.8 * const.milli
f0  = 1.35 * const.giga
f1  = 1.25 * const.giga
f2  = 1.45 * const.giga

Z1  = np.sqrt( 2) * Z0
Z3  = Z0 / np.sqrt( 2)

print( 'Z1 = {}; Z3 = {}' .format( Z1, Z3))

def find_w( Z0, b, t, mu, eps):
    m  = 6 * ( b - t) / ( 3 * b - t)
    B  = np.exp( Z0 * np.sqrt( eps) / 30)
    W  = 8 * ( ( b - t) / const.pi) * ( np.sqrt( B + 0.568) / ( B - 1))
    dw = t * \
         ( 1 - \
           ( np.log( ( ( t / ( 2 * b - t)) ** 2) + \
                     ( ( ( 0.0796 * t) / ( W - (0.26*t))) ** m) / 2))) / const.pi
    return W - dw

print 'wt = {}' .format( find_w( 30, \
                                 2*const.milli, \
                                 0.01*const.milli, \
                                 1, \
                                 10.20) / const.milli)

print 'w0 = {}' .format( find_w( Z0, b, t, 1, 1) / const.milli)
print 'w1 = {}' .format( find_w( Z1, b, t, 1, 1) / const.milli)
print 'w3 = {}' .format( find_w( Z3, b, t, 1, 1) / const.milli)

def coupling(f, f0):
    theta  = ( const.pi / 2) * f / f0
    c = np.cos( theta)
    s = np.sin( theta)
    t = np.tan( theta)

#    g1 = ( 1 + ( t ** 2)) / ( t ** 2) + 3 + ( 2 * np.sqrt( 2))
    be = ( np.sqrt( 2) * t - ( ( 2 + np.sqrt( 2)) * ( t ** 3))) / \
         ( 2 * ( t ** 4) - ( ( 2 * np.sqrt( 2) - 1) * ( t ** 2)) + 1)
    ge = ( 1 + ( t ** 2)) / \
         ( 2 * ( t ** 4) - ( 2 * np.sqrt( 2) - 1) * ( t ** 2) + 1)
    R22 = c * ( 3 + ( 2 * ge)) - ( s * be * np.sqrt( 2))
    X22 = ( s * ge * np.sqrt( 2)) + \
          ( s * 2 * np.sqrt( 2)) + \
          ( 2 * be * c)
    S12 = 2 / ( R22 + ( 1j * X22))
    return 20 * np.log10( 1 / np.absolute( S12))

freq = np.linspace( f1, f2)
line, = plt.plot( freq, coupling( freq, f0))
plt.xlabel( 'czestotliwosc [Hz]')
plt.ylabel( 'izolacja I [dB]')
plt.show()
