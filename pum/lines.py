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

def cylindrical_flat_coupled( s, d, h, mu, epsilon):
    x = d / h
    y = s / h
    a = 1 + np.exp( 16 * x - 18.272)
    b = np.sqrt( 5.905 - ( x**4))
    c = ( -0.8107 * ( y**3)) + \
        ( 1.3401 * ( y**2)) - \
        ( 0.6929 * y) + \
        ( 1.0892) + \
        ( 0.014002 / y) - \
        ( 0.000636 / ( y**2))
    d = 0.11 - ( 0.83 * y) + ( 1.64 * ( y**2)) - ( y**3)
    e = -0.15 * np.exp( -13 * x)
    g = 2.23 * np.exp( ( -7.01 * y) + \
                       ( 10.24 * ( y**2)) - \
                       ( 27.58 * ( y**3)))
    k = 1 + ( 0.01 * ( ( -0.0726) - \
                       ( 0.2145 / y) + \
                       ( 0.222573 / ( y ** 2)) - \
                       ( 0.012823 / ( y ** 3))))
    l = 0.01 * ( ( -0.26) + \
                 ( 0.6866 / y) + \
                 ( 0.0831 / ( y ** 2)) - \
                 ( 0.0076 / ( y ** 3)))
    m = ( -0.1098) + \
        ( 1.2138 * x) - \
        ( 2.2535 * ( x ** 2)) + \
        ( 1.1313 * ( x ** 3))
    n = ( -0.019) - \
        ( 0.016   / y) + \
        ( 0.0362  / ( y ** 2)) - \
        ( 0.00234 / ( y ** 3))
    f1 = x * a / b
    if y < 0.9:
        f2 = c - ( x * d) + ( e * g)
    else:
        f2 = 1 + ( 0.004 * np.exp( 0.9 - y))
    f3 = np.tanh( const.pi * ( x + y) / 2)
    if y < 0.9:
        f4 = k - ( x * l) + ( m * n)
    else:
        f4 = 1
    Z0e = 59.952 * np.log( 0.523962 / ( f1 * f2 * f3))
    Z0o = 59.952 * np.log( ( 0.523962 * f3) / ( f1 * f4))
    return ( Z0e, Z0o)
