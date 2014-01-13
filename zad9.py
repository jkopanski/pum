import scipy.constants as const
import scipy.optimize  as opt
import numpy           as np
import pum.algorithms  as alg
from pum.lines import *
from pum.net import *

L   = 10.0
Z01 = 50.0
Z02 = 60.0

r = Z01 / Z02
if r <= 1:
    r = 1 / r
print 'r = {0}' .format( r)

Lmin = 10 * np.log( np.sqrt( r) + np.sqrt( r - 1))
print 'Lmin = {0}' .format( Lmin)

N  = 10 ** ( L / 10)
print 'N = {0}' .format( N)

R3 = ( 2 * np.sqrt( N * Z01 * Z02)) / ( N - 1)
R1 = Z01 * ( N + 1) / ( N - 1) - R3
R2 = Z02 * ( N + 1) / ( N - 1) - R3

print 'R1 = {0}; R2 = {1}; R3 = {2}' .format( R1, R2, R3)

Ra = ( ( N - 1) * np.sqrt( Z01 * Z02)) / ( 2 * np.sqrt(N))
Rb = ( Z01 * Ra * ( N - 1)) / ( ( Ra * ( N + 1)) - ( Z01 * ( N - 1)))
Rc = ( Z02 * Ra * ( N - 1)) / ( ( Ra * ( N + 1)) - ( Z02 * ( N - 1)))
print 'Ra = {0}; Rb = {1}; Rc = {2}' .format( Ra, Rb, Rc)
