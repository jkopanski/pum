import scipy.constants as const
import numpy as np

def coax_z(a, b, mu, epsilon): # impedance of coaxial line
    return np.sqrt(const.mu_0 * mu / ( const.epsilon_0 * epsilon)) * np.log( a / b) / ( 2 * const.pi)

def skew_coax_z( c, a, b, mu, epsilon):
    x = ( b + ( a * a - ( 4 * c * c)) / b) / (2 * a)
    return 59.952 * np.sqrt( mu / epsilon) * np.log( x + np.sqrt( x * x - 1))

def cylindrical_flat( d, b, mu, epsilon): # impedance of cylindrical flat line
    R = const.pi * d / (4 * b)
    x = 1 + 2 * np.sinh( R) ** 2
    y = 1 - 2 * np.sin( R) ** 2
    return 59.952 * np.sqrt( mu / epsilon) * ( np.log( np.sqrt( x) + np.sqrt( y) / np.sqrt( x - y)) - R**4 / 30 + 0.014 * ( R**8))
