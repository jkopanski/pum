import scipy.constants as const
import scipy.integrate as integ
import numpy as np
import matplotlib.pyplot as plt
import shapely.geometry as shp

class pt:
    """Single point of finite difference method net"""
    def __init__( self, u = 0.0, p = 1.0, q = 1.0, r = 1.0):
        self.__dict__['u'] = float( u)
        self.p = float( p)
        self.q = float( q)
        self.r = float( r)

class metal( pt):
    """Metal point.
    
    Assign new potential only upon creation"""
    def __setattr__( self, name, val):
        if name == 'u':
            pass
        else:
            self.__dict__[name] = val

class imaginary( pt):
    """Imaginary point.

    Used when dividing structure"""
    def __init__( self, pt):
        self.ref = pt

    def __getattr__( self, name):
        if name == 'u':
            return self.ref.u

class struct:
    """Microwave structure to analyze with finite difference method"""
    def __init__( self, a, b, I, J, mu=1, eps=1):
        self.a = a
        self.b = b
        self.I = I
        self.J = J
        self.mu = mu
        self.eps = eps
        self.xstep = a / I
        self.ystep = b / J
        self.st = [[pt(0.0) for x in xrange(self.I+1)] \
                   for z in range(self.J+1)]
        self.bl = []

    def add_rect( self, a, b, x, y, u):
        xs = a / 2
        ys = b / 2
        self.bl.append( ( shp.Polygon( [ (x - xs, y + ys), \
                                         (x + xs, y + ys), \
                                         (x + xs, y - ys), \
                                         (x - xs, y - ys)]), u))

    def add_plane( self, sx, sy, ex, ey, u):
        self.bl.append( ( shp.LineString( [(sx, sy), (ex, ey)]), u))

    def init( self):
        for j in range( self.J + 1):
            for i in range( self.I + 1):
                # is it metal?
                for item in self.bl:
                    #if item[0].contains( shp.Point( i * self.xstep, j * self.ystep)):
                    if not item[0].disjoint( shp.Point( i * self.xstep, j * self.ystep)):
                        self.st[j][i] = metal( item[1])
                        continue
                if not isinstance( self.st[j][i], metal):
                    # if not is is standard point
                    # but look for metal vicinity
                    for item in self.bl:
                        p = self.ystep / self.xstep
                        q = 1.0
                        r = p
                        step = self.xstep
                        #if item[0].contains( shp.Point( (i - 1) * self.xstep, (j + 0) * self.ystep)):
                        if not item[0].disjoint( shp.Point( (i - 1) * self.xstep, (j + 0) * self.ystep)):
                            # calculate new p, q, r values
                            line = shp.LineString( [( i * self.xstep, j * self.ystep), ( (i - 1) * self.xstep, (j + 0) * self.ystep)])
                            inter = item[0].intersection( line)
                            step = shp.Point( i * self.xstep, j * self.ystep).distance( inter)
                            p = p * self.xstep / step
                            q = q * self.xstep / step
                            r = r * self.xstep / step
                        #if item[0].contains( shp.Point( (i + 1) * self.xstep, (j + 0) * self.ystep)):
                        if not item[0].disjoint( shp.Point( (i + 1) * self.xstep, (j + 0) * self.ystep)):
                            # calculate new q value
                            line = shp.LineString( [( i * self.xstep, j * self.ystep), ( (i + 1) * self.xstep, (j + 0) * self.ystep)])
                            inter = item[0].intersection( line)
                            #q = ( inter.bounds[0] - ( i * self.xstep)) / step
                            q = shp.Point( i * self.xstep, j * self.ystep).distance( inter) / step
                        #if item[0].contains( shp.Point( (i + 0) * self.xstep, (j + 1) * self.ystep)):
                        if not item[0].disjoint( shp.Point( (i + 0) * self.xstep, (j + 1) * self.ystep)):
                            # calculate new p value
                            line = shp.LineString( [( i * self.xstep, j * self.ystep), ( (i + 0) * self.xstep, (j + 1) * self.ystep)])
                            inter = item[0].intersection( line)
                            #p = ( inter.bounds[1] - ( j * self.ystep)) / step
                            p = shp.Point( i * self.xstep, j * self.ystep).distance( inter) / step
                        #if item[0].contains( shp.Point( (i + 0) * self.xstep, (j - 1) * self.ystep)):
                        if not item[0].disjoint( shp.Point( (i + 0) * self.xstep, (j - 1) * self.ystep)):
                            # calculate new r value
                            line = shp.LineString( [( i * self.xstep, j * self.ystep), ( (i + 0) * self.xstep, (j - 1) * self.ystep)])
                            inter = item[0].intersection( line)
                            #r = ( ( j * ystep) - inter.bounds.maxy) / step
                            r = shp.Point( i * self.xstep, j * self.ystep).distance( inter) / step
                    #add point
                    self.st[j][i] = pt( 0.0, p, q, r)

    def calc(self, i, j):
        u =  self.st[j + 0][i - 1].u / \
             ( 1.0 + self.st[j][i].q)
        u += self.st[j - 0][i + 1].u / \
             ( self.st[j][i].q * ( 1.0 + self.st[j][i].q))
        u += self.st[j + 1][i - 0].u / \
             ( self.st[j][i].p * (self.st[j][i].p + self.st[j][i].r))
        u += self.st[j - 1][i + 0].u / \
             ( self.st[j][i].r * (self.st[j][i].p + self.st[j][i].r))   
        u *= ( self.st[j][i].p * self.st[j][i].q * self.st[j][i].r)
        u /= ( self.st[j][i].p * self.st[j][i].r + self.st[j][i].q)
        return u

    def liebmann(self, tol=1e-8):
        k = 0
        while True:
            Rm = 0.0
            k += 1
            for j in range( self.J + 1):
                for i in range( self.I + 1):
                    #sprawdzic jaki to punkt
                    if not isinstance( self.st[j][i], metal):
                        v = self.st[j][i].u
                        self.st[j][i].u = self.calc( i, j)
                        R = np.abs(self.st[j][i].u - v)
                        if Rm < R:
                            Rm = R
            if Rm < tol:
                break
        return k

    def liebmann_quater(self, tol=1e-8):
        k = 0
        for j in range( self.J/2 + 1):
            self.st[j][self.I/2+1] = imaginary( self.st[j][self.I/2-1])
        for i in range( self.I/2 + 1):
            self.st[self.J/2+1][i] = imaginary( self.st[self.J/2-1][i])
        while True:
            Rm = 0.0
            k += 1
            for j in range( self.J/2+1):
                for i in range( self.I/2+1):
                    #sprawdzic jaki to punkt
                    if not ( isinstance( self.st[j][i], metal) and \
                             isinstance( self.st[j][i], imaginary)):
                        v = self.st[j][i].u
                        self.st[j][i].u = self.calc( i, j)
                        R = np.abs(self.st[j][i].u - v)
                        if Rm < R:
                            Rm = R
            if Rm < tol:
                break
        for j in range( self.J + 1):
            for i in range( self.I + 1):
                if j < self.J / 2 + 1 and i > self.I / 2:
                    self.st[j][i] = self.st[j][self.I - i]
                elif j > self.J / 2:
                    self.st[j][i] = self.st[self.J - j][i]
        return k

    def impedance( self):
        Eyt = [float for x in xrange( self.I+1)]
        Eyd = [float for x in xrange( self.I+1)]
        Exl = [float for x in xrange( self.J+1)]
        Exr = [float for x in xrange( self.J+1)]

        for j in range( self.J + 1):
            # left boundary
            d1 =      ( self.st[j][0].q * self.xstep)
            d2 = d1 + ( self.st[j][1].q * self.xstep)
            e1st = ( self.st[j][1].u - self.st[j][0].u) / d1
            e2nd = ( self.st[j][2].u - self.st[j][0].u) / d2
            Exl[j] = e1st + ( ( e1st - e2nd) / ( ( d2 / d1) - 1))
            # right boundary
            d1 =      ( self.st[j][self.I-1].q * self.xstep)
            d2 = d1 + ( self.st[j][self.I-2].q * self.xstep)
            e1st = ( self.st[j][self.I-1].u - self.st[j][self.I-0].u) / d1
            e2nd = ( self.st[j][self.I-2].u - self.st[j][self.I-0].u) / d2
            Exr[j] = e1st + ( ( e1st - e2nd) / ( ( d2 / d1) - 1))

        for i in range( self.I + 1):
            # lower boundry
            d1 =      ( self.st[0][i].p * self.xstep)
            d2 = d1 + ( self.st[1][i].p * self.xstep)
            e1st = ( self.st[1][i].u - self.st[0][i].u) / d1
            e2nd = ( self.st[2][i].u - self.st[0][i].u) / d2
            Eyd[i] = e1st + ( ( e1st - e2nd) / ( ( d2 / d1) - 1))
            # upper boundry
            d1 =      ( self.st[self.J-0][i].r * self.xstep)
            d2 = d1 + ( self.st[self.J-1][i].r * self.xstep)
            e1st = ( self.st[self.J-1][i].u - self.st[self.J-0][i].u) / d1
            e2nd = ( self.st[self.J-2][i].u - self.st[self.J-0][i].u) / d2
            Eyt[i] = e1st + ( ( e1st - e2nd) / ( ( d2 / d1) - 1))
            
        E = integ.simps( Exl, None, self.ystep) + \
            integ.simps( Exr, None, self.ystep) + \
            integ.simps( Eyt, None, self.xstep) + \
            integ.simps( Eyt, None, self.xstep)
        # print 'Exl = {}; Exr = {}; Eyt = {}; Eyb = {}' \
        #     .format( integ.simps( Exl, None, self.ystep), \
        #              integ.simps( Exr, None, self.ystep), \
        #              integ.simps( Eyt, None, self.xstep), \
        #              integ.simps( Eyt, None, self.xstep))
        # print 'E = {}' .format( E)
        return np.sqrt( ( const.mu_0 * self.mu) / \
                        ( const.epsilon_0 * self.eps)) * \
            ( 1 / E)

    def plot( self):
        z = [[float for x in xrange(self.I+1)] for z in range(self.J+1)]
        for j in range( self.J + 1):
            for i in range( self.I + 1):
                z[j][i] = self.st[j][i].u
        cmap = plt.get_cmap( 'jet')
        plt.imshow( z, origin='lower', interpolation='nearest', cmap=cmap)
        plt.colorbar()
        plt.show()
