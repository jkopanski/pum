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
plt.show()
