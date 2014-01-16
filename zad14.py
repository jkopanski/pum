import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
import matplotlib.pyplot as plt
from pum.lines import *
from pum.net import *

Z01 = 35.0
Z0  = 50.0
f0  = 1.35 * const.giga
f1  = 1.25 * const.giga
f2  = 1.45 * const.giga

q   = Z01 / Z0
print( 'q = {}' .format( q))

w  = ( 2 * (f2 - f1)) / (f1 + f2)
u0 = np.sin( w * const.pi / 4)
C  = ( ( ( 2 * q) - 1) * ( u0 ** 2)) / \
     ( 2 * ( 2 - ( u0 ** 2)))
v1sq = np.sqrt( ( C ** 2) + ( 2 * q)) + C
v2 = 2 * q / v1sq

Z1 = Z0 * np.sqrt( v1sq)
Z2 = Z1 * v2

print( 'Z1 = {}; Z2 = {}' .format( Z1, Z2))

theta3 = ( const.pi / 2) * ( 1 - ( w / ( 2 * np.sqrt( 2))))
R2 = ( 2 * Z1 * Z2) / \
     np.sqrt( ( Z1 + Z2) * \
              ( Z2 - Z1 * ( ( 1 / np.tan( theta3)) ** 2)))
R1 = ( 2 * R2 * (Z1 + Z2)) / \
     ( R2 * ( Z1 + Z2) / Z0 - ( 2 * Z2))

print( 'R1 = {}; R2 = {}' .format( R1, R2))

def swr(f, f0, Z0, Z1, Z2, R1, R2):
    z1 = Z1 / Z0
    z2 = Z2 / Z0
    y1 = 1 / z1
    y2 = 1 / z2
    g1 = Z0 / R1
    g2 = Z0 / R2
    theta  = ( const.pi / 2) * f / f0

    A = ( y1 + y2) * ( 1 - ( 2 * g1)) - ( 2 * g2 * y1)
    B = ( 2 * g2 - ( 4 * g1 * g2) - ( y1 ** 2)) * \
        np.tan( theta) + \
        ( y1 * y2 * ( 1 / np.tan( theta)))
    C = ( y1 + y1) * ( 1 + ( 2 * g1)) + ( 2 * g2 * y1)
    D = ( 2 * g2 + ( 4 * g1 * g2) + ( y1 ** 2)) * \
        np.tan( theta) - \
        ( y1 * y2 * ( 1 / np.tan( theta)))
    Spm = ( A + ( 1j * B)) / ( C + ( 1j * D))

    RD = ( ( 2 * q + 1) * ( np.cos( theta) ** 2)) - \
         ( ( ( z2 * y1) + ( 2 * q * z1 * y2)) * ( np.sin( theta) ** 2))
    RN = ( ( 2 * q - 1) * ( np.cos( theta) ** 2)) + \
         ( ( ( z2 * y1) - ( 2 * q * z1 * y2)) * ( np.sin( theta) ** 2))
    XD = ( z1 + z2 + ( 2 * q * ( y1 + y2))) * \
         np.cos( theta) * np.sin( theta)
    XN = ( z1 + z2 - ( 2 * q * ( y1 + y2))) * \
         np.cos( theta) * np.sin( theta)
    S11 = ( -RN + ( 1j * XN)) / ( RD + ( 1j * XD))
    return ( 1.0 + np.absolute( S11)) / ( 1.0 - np.absolute( S11))

def isolation(f, f0, Z0, Z1, Z2, R1, R2):
    z1 = Z1 / Z0
    z2 = Z2 / Z0
    y1 = 1 / z1
    y2 = 1 / z2
    g1 = Z0 / R1
    g2 = Z0 / R2
    theta  = ( const.pi / 2) * f / f0

    A = ( y1 + y2) * ( 1 - ( 2 * g1)) - ( 2 * g2 * y1)
    B = ( 2 * g2 - ( 4 * g1 * g2) - ( y1 ** 2)) * \
        np.tan( theta) + \
        ( y1 * y2 * ( 1 / np.tan( theta)))
    C = ( y1 + y1) * ( 1 + ( 2 * g1)) + ( 2 * g2 * y1)
    D = ( 2 * g2 + ( 4 * g1 * g2) + ( y1 ** 2)) * \
        np.tan( theta) - \
        ( y1 * y2 * ( 1 / np.tan( theta)))
    Spm = ( A + ( 1j * B)) / ( C + ( 1j * D))

    RD = ( ( 2 * q + 1) * ( np.cos( theta) ** 2)) - \
         ( ( ( z2 * y1) + ( 2 * q * z1 * y2)) * ( np.sin( theta) ** 2))
    RN = ( ( 2 * q - 1) * ( np.cos( theta) ** 2)) + \
         ( ( ( z2 * y1) - ( 2 * q * z1 * y2)) * ( np.sin( theta) ** 2))
    XD = ( z1 + z2 + ( 2 * q * ( y1 + y2))) * \
         np.cos( theta) * np.sin( theta)
    XN = ( z1 + z2 - ( 2 * q * ( y1 + y2))) * \
         np.cos( theta) * np.sin( theta)
    S23 = ( ( RN + ( 1j * XN)) / ( RD + ( 1j * XD)) - Spm) / 2
    return 20 * np.log10( 1 / np.absolute( S23))

freq = np.linspace( f1, f2)
line, = plt.plot( freq, swr( freq, f0, Z0, Z1, Z2, R1, R2))
line2, = plt.plot( freq, isolation( freq, f0, Z0, Z1, Z2, R1, R2))
plt.show()
