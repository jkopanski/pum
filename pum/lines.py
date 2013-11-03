from scipy import constants
import numpy as np

def coax_z(a, b, mu, epsilon): # impedance of coaxial line
    return np.sqrt(constants.mu_0 * mu / ( constants.epsilon_0 * epsilon)) * np.log( a / b) / ( 2 * constants.pi)

def skew_coax_z( c, a, b, mu, epsilon):
    x = ( b + ( a * a - ( 4 * c * c)) / b) / (2 * a)
    return 59.952 * np.sqrt( mu / epsilon) * np.log( x + np.sqrt( x * x - 1))
