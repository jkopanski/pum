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

def microstrip(w, t, h, f, mu, epsilon):
    u  = w / h
    if t != 0:
        du = ( t / ( 2 * const.pi * h)) * np.log( 1 + ( ( 4 * np.exp(1) * h) / ( t * ( ( 1 / np.tanh( np.sqrt( 6.517 * u)))**2)))) * ( 1 + ( 1 / np.cosh( np.sqrt( epsilon - 1))))
        u += du
    a  = 1 + ( ( 1 / 49) * np.log( ( ( u**4) + ( ( u / 52)**2)) / ( ( u**4) + 0.432))) + ( ( 1 / 18.7) * np.log( 1 + ( ( u / 18.1)**3)))
    b  = 0.564 * ( ( ( epsilon - 0.9) / ( epsilon + 3))**0.053)
    c  = 0.33622 * ( 1 - np.exp( -0.03442 * epsilon))
    d  = 3.751 - ( 2.75 * np.exp( - ( ( epsilon / 15.916)**8)))
    fn = f * h * ( 10**-7)
    g  = 1 - np.exp( - ( ( fn / 3.87)**4.97))
    k  = 0.363 * g * np.exp( -4.6 * u)
    m  = 0.525 / ( ( 1 + ( 0.157 * fn))**20)
    n  = 0.27488 + ( 0.6315 * u) + ( m * u) - ( 0.065683 * np.exp( -8.7513 * u))
    p  = n * c * ( ( ( 1.844 * fn) + ( k * d * fn))**1.5763)
    eps_0 = ( ( epsilon + 1) / 2) + ( ( ( epsilon - 1) / 2) * ( ( 1 + ( 10 / u))**( -a * b)))
    eps_eff = ( eps_0 + ( epsilon * p)) / ( 1 + p)
    f = 6 + ( ( 2 * const.pi - 6) * np.exp( - ( ( 30.666 / u)**0.7528)))
    return ( 60 / np.sqrt( eps_eff)) * np.log( ( f / u) + np.sqrt( 1 + ( 2 / u)**2))
