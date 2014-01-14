import numpy as np

def newton_raphson( func, a, b, args=(), tol=1e-8):
    """
    Znajduje miejsce zerowe funkcji korzystajac z metody Newtona-Raphsona

    Parameters
    ----------
    func : function
        Funkcja ktorej miejsce zerowe jest poszukiwane
    a : float
        Poczatek zakresu na ktorym szukane jest miejsce zerowe
    b : float
        Koniec zakresu na ktorym szukane jest miejsce zerowe
    args : tuple, optional
        Dodatkowe parametry przekazywane do funkcji
    tol : float, optional
        Dopuszczalny blad znalezionego moejsca zerowego

    Returns
    ----------
    zero : float
        Szacowana lokacja gdzie wartosc funkcji wynosi 0.
    """
    ###
    # Szukanie punktu startowego
    ###
    x0 = a * 1.0 # mnozenie przez 1.0 sprawia, ze zmienna staje sie float
    delta_x = ( b - a) / 25
    myargs = (x0,) + args
    fval_min = np.abs( func( *myargs))
    zero = x0
    for x0 in np.arange( a, b, delta_x):
        myargs = (x0,) + args
        fval = np.abs( func( *myargs))
        if fval < fval_min:
            fval_min = fval
            zero = x0
            
    ###
    # Wlasciwa metoda newtona
    ###
    delta_x = 1e-8
    error = 1
    while error > tol:
        x = zero
        fd1 = func( *(x - delta_x,) + args)
        fd2 = func( *(x + delta_x,) + args)
        fderiv = ( fd2 -fd1) / ( 2 * delta_x)
        if fderiv == 0:
            msg = "derivative was zero."
            warnings.warn(msg, RuntimeWarning)
        fval = func( *(x,) + args)
        zero = x - ( fval / fderiv)
        error = zero - x

    return zero

def n_fun( q, tol=1e-19):
    """
    Zwraca wartosc N potrzebna do wyznaczenia modulu

    Parameters
    ----------
    q : float
        stala modularna
    tol : float, optional
        dokladnosc

    Returns
    ----------
    n : float
        Wartosc parametru N.
    """
    n = 0.0
    old = -1.0
    i = 1
    while( n - old > tol):
        old = n
        n = old + q ** ( i * ( i - 1))
        i = i + 1
    return n

def d_fun( q, tol=1e-19):
    """
    Zwraca wartosc D potrzebna do wyznaczenia modulu

    Parameters
    ----------
    q : float
        stala modularna
    tol : float, optional
        dokladnosc

    Returns
    ----------
    d : float
        Wartosc parametru D.
    """
    d = 0.5
    old = 0.0
    i = 1
    while( d - old > tol):
        old = d
        d = old + q ** ( i * i)
        i = i + 1
    return d

def k_int( k, tol=1e-8):
    old = 0
    k0 = k
    k0p = np.sqrt( 1 - ( k ** 2))
    k_int = ( 1 + k0) / ( 1 + k0p)
    i = 1
    while( k_int - old > tol):
        old = k_int
        k0  = ( 2 / ( 1 + k0 )) * np.sqrt( k0)
        k0p = ( 2 / ( 1 + k0p)) * np.sqrt( k0p)
        k_int = k_int * ( 1 + k0) / ( 1 + k0p)
        i += 1
    return k_int
