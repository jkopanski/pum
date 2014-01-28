import numpy as np
import matplotlib.pyplot as plt

class point:
    """Klasa reprezentujaca podstawowy punkt siatki"""
    def __init__(self, u=0.0, p=1.0, q=1.0, r=1.0):
        self.u = float( u)
        self.p = float( p)
        self.q = float( q)
        self.r = float( r)

    def update(self, u):
        self.u = u

class met_point( point):
    def __init__(self, u=0.0, p=1.0, q=1.0, r=1.0):
        self.__dict__['u'] = float( u)
        self.p = float( p)
        self.q = float( q)
        self.r = float( r)

    def update(self, u):
        pass

    def __setattr__( self, name, val):
        if name == 'u':
            pass
        else:
            self.__dict__[name] = val
        
class imag_point( point):
    def __init__(self, pt):
        self.ref = pt

    def __getattr__( self, name):
        if name == 'u':
            return self.ref.u

    def update(self, u):
        self.u = self.ref.u


class net:
    """Klasa reprezentujaca siatke do obliczen metoda roznic skonczonych"""
    # def __init__(self, a, b, t, w, I, J):
    #     self.a = a
    #     self.b = b
    #     self.t = t
    #     self.w = w
    #     self.I = I
    #     self.J = J
    #     self.x = I + 1
    #     self.y = J + 1
    #     self.net = [point for x in range( self.x * self.y)]
    #     for j in range( self.x + 1):
    #         for i in range( self.y + 1):
    #             if j == 0:
    #                 self.net[j * (self.x + 1) + i] = met_point(0.0)
    #             elif i == 0:
    #                 self.net[j * (self.x + 1) + i] = met_point(0.0)
    #             else:
    #                 self.net[j * (self.x + 1) + i] = point(0.0)
    def __init__(self, a, b, I, J, Is, Js, q=1, p=1, r=1):
        self.a  = float(a)
        self.b  = float(b)
        self.I  = I
        self.J  = J
        self.Is = Is
        self.Js = Js
        self.q  = float(q)
        self.p  = float(p)
        self.r  = float(r)
#        self.net = [point(0.0) for x in range((I) * (J))]
        self.net = [[point(0.0) for x in xrange(self.I)] for z in range(self.J)]

        print 'I = {}; J = {}' .format( self.I, self.J)

    def init(self):
        for j in range( 1, self.J-1):
            for i in range( 1, self.I-1):
                p  = (self.b / self.a) * (float(self.I-1)/float(self.J-1)) * (float(self.J - 1 - j)/float(i))
                q  = float( self.I - 1 - i) / float(i)
                r  = (self.b / self.a) * (float(self.I-1)/float(self.J-1)) * (float(j)/float(i))
                u  = self.net[j][0].u / ( 1 + q)
                u += self.net[j][self.I - 2].u / ( q * (1 + q))
                u += self.net[self.J-1][i].u / ( p * (p + r))
                u += self.net[0][i].u / ( r * (p + r))
                u *= ( p * q * r)
                u /= ( p * r + q)
                self.net[j][i].u = u

    def calc(self, i, j):
        u =  self.net[j + 0][i - 1].u / (            1.0 + self.q)
        u += self.net[j - 0][i + 1].u / ( self.q * (   1.0 + self.q))
        u += self.net[j + 1][i - 0].u / ( self.p * (self.p + self.r))
        u += self.net[j - 1][i + 0].u / ( self.r * (self.p + self.r))                
        # u =  self.net[(j + 0)]*self.I + i - 1].u / (            1.0 + self.q)
        # u += self.net[(j - 0)*self.I + i + 1].u / ( self.q * (   1.0 + self.q))
        # u += self.net[(j + 1)*self.I + i - 0].u / ( self.p * (self.p + self.r))
        # u += self.net[(j - 1)*self.I + i + 0].u / ( self.r * (self.p + self.r))                
        u *= ( self.p * self.q * self.r)
        u /= ( self.p * self.r + self.q)
        return u

    def liebmann(self, tol=1e-8):
        k = 0
        while True:
            Rm = 0.0
            k += 1
            for j in range( 1, self.J - 1):
                for i in range( 1, self.I - 1):
                    #sprawdzic jaki to punkt
                    if type( self.net[j][i]) != met_point and \
                       type( self.net[j][i]) != imag_point:
                        v = self.net[j][i].u
                        self.net[j][i].u = self.calc( i, j)
                        R = np.abs(self.net[j][i].u - v)
                        if Rm < R:
                            Rm = R
            if Rm < tol:
                break
        return k

    def plot(self):
        k = 0
        z = [[0.0 for x in xrange(self.I)] for z in range(self.J)]
        for j in range( self.J):
            for i in range( self.I):
                z[j][i] = self.net[j][i].u
                k += 1
        x = np.mgrid[0:self.I]
        y = np.mgrid[0:self.J]
        cmap = plt.get_cmap( 'gist_stern')
        plt.imshow( z, origin='lower', interpolation='nearest', cmap=cmap)
        plt.colorbar()
        plt.show()
