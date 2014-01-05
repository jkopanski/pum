class net:
    """Klasa reprezentujaca siatke do obliczen metoda roznic skonczonych"""
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
        self.net = [float(0) for x in range((I) * (J))]
        print '{}' .format(self.net)
        
    def init(self):
        for j in range( 1, self.J-1):
            for i in range( 1, self.I-1):
                p  = (self.b / self.a) * (float(self.I-1)/float(self.J-1)) * (float(self.J - 1 - j)/float(i))
                q  = float( self.I - 1 - i) / float(i)
                r  = (self.b / self.a) * (float(self.I-1)/float(self.J-1)) * (float(j)/float(i))
                u  = self.net[j*self.I + 0] / ( 1 + q)
                u += self.net[j*self.I + self.I-1] / ( q * (1 + q))
                u += self.net[(self.J-1)*self.I + i] / ( p * (p + r))
                u += self.net[0*self.I + i] / ( r * (p + r))
                u *= ( p * q * r)
                u /= ( p * r + q)
                self.net[j*self.I + i] = u

    def calc(self, i, j):
        u =  self.net[(j + 0)*self.I + i - 1] / (                1 + self.q)
        u += self.net[(j - 0)*self.I + i + 1] / ( self.q * (     1 + self.q))
        u += self.net[(j + 1)*self.I + i - 0] / ( self.p * (self.p + self.r))
        u += self.net[(j - 1)*self.I + i + 0] / ( self.r * (self.p + self.r))                
        # stare
        # u =  self.net[j*self.I + (i-1)] / ( 1 + self.q)
        # u += self.net[j*self.I + (i+1)] / ( self.q * (     1 + self.q))
        # u += self.net[(j+1)*self.I + i] / ( self.p * (self.p + self.r))
        # u += self.net[(j-1)*self.I + i] / ( self.r * (self.p + self.r))
        u *= ( self.p * self.q * self.r)
        u /= ( self.p * self.r + self.q)
        return u

    def liebmann(self, tol=1e-8):
        Rm = 0
        k = 0
        while True:
            k += 1
            print 'k = {}'  .format( k)
            print 'Rm = {}' .format( Rm)
            print '{}'      .format( self.net)
            for j in range( 1, self.J - 1):
                for i in range( 1, self.I - 1):
                    v = self.net[j*self.I + i]
                    self.net[j*self.I + i] = self.calc( i, j)
                    R = abs( self.net[j*self.I + i] - v)
                    if Rm < R:
                        Rm = R
            if k == 16:
                break
        return k
