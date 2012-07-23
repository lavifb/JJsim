import math, random
from DiffSolver import euler, matmult, matInv

class JJ:
    """ A Josephson Junction"""
    def __init__(self, b_c = 1.0, p = 0.0, v = 0.0, dt = .01):
        """ Initiates the junction.

            b_c: Stewart-McCumber damping constant
            p: the initial phase difference across the junction
            v: average voltage across the junction
            dt: timestep """
        self.phase = float(p)
        self.volt = float(v)
        self.b = float(b_c)
        # b>1 is underdmaped. b<1 is overdamped
        self.dt = float(dt)
        self.i = 0

    def getInfo(self):
        """ Returns info about junction."""
        out = 'b = {0:n}'.format(self.b)
        return out

    def getType(self):
        """ Returns type of junction."""
        return 'simple'

    def applyI(self, i = 0.0, T=1000):
        """ Applies bias current to the junction for duration T.

            i: applied bias current
            dt: timestep for numerics
            T: duration of counting

            sets the phase to the most recent phase
            sets voltage to average over last three fifths 
            of T returns the average voltage."""
        self.i = i
        iT = int(T/self.dt) # renormalized time to integer values
        p = self.phase
        v = self.volt
        sumv = 0.0
        vtot = 0
        for t in xrange(iT):
            pv = euler([p, v], self.dx, self.dt)
            p, v = pv
            if t > 3*iT/5:
                sumv+=v
                vtot+=1
        self.setP(p)
        self.volt = v
        return sumv/vtot

    def dx(self, x, t):
        """ Returns derivative for diff eq."""
        i_ = self.i
        b_ = self.b
        dp = x[1]
        dv = (i_ - x[1] - math.sin(x[0]))/b_
        return [dp, dv]

    def setP(self, p):
        """ Sets the phase, taking it mod 2*pi"""
        self.phase = math.fmod(p, 2*math.pi)
        pass

    def getTemp(self):
        """ Returns normalized temperature at which the junction is operating."""
        return 0

    def getPhaseVolt(self, i = 0.0):
        """ Applies current and returns phase and voltage at point."""
        self.i = i
        pv = euler([self.phase, self.volt], self.dx, self.dt)
        self.phase, self.volt = pv
        return pv

class JJn(JJ):
    """A Josephson Junction with thermal noise"""
    def __init__(self, b_c = 1.0, p = 0.0, v = 0.0, dt = .01, temp = 0.0):
        """ Initiates the junction.

            b_c: Stewart-McCumber damping constant
            p: the initial phase difference across the junction
            v: average voltage across the junction
            dt: timestep 
            temp: unit-less temperature """

        JJ.__init__(self, b_c, p, v, dt)
        self.setTemp(temp)

    def getInfo(self):
        """ Returns info about junction."""
        out = 'b = {0:n}, temp = {1:n}'.format(self.b, self.getTemp())
        return out

    def getType(self):
        """ Returns type of junction"""
        return 'noisy'

    def dx(self, x, t):
        """ Returns derivative for diff eq."""
        i_ = self.i
        b_ = self.b
        i_n = random.gauss(0, self.sig)
        dp = x[1]
        dv = (i_ + i_n - x[1] - math.sin(x[0]))/b_
        return [dp, dv]

    def getTemp(self):
        """ Returns normalized temperature at which the junction is operating."""
        return self.sig*self.sig*.5*self.dt

    def setTemp(self, temp):
        """ Sets self.sig to the appropriate value for the temperature."""
        self.sig = math.sqrt((temp + temp)/self.dt)

class JJFreq(JJ):
    """ A Josephson Junction with frequency dependent circuit elements."""
    def __init__(self, b_c = 1.0, p = 0.0, v = 0.0, vc = 0.0, Q1 = 1.0, rho = 1.0, dt = .01):
        """ Initiates the junction.

            b_c: Stewart-McCumber damping constant
            p: the initial phase difference across the junction
            v: average voltage across the junction
            dt: timestep 
            vc: average voltage across the capacitor
            Q1: quality factor of the frequency dependent branch of the junction
            rho: quotient of the different junction time constants. """

        JJ.__init__(self, b_c, p, v, dt)
        self.v_c = float(vc)
        self.d = math.sqrt(self.b)/Q1 - 1
        self.e = rho/self.b

    def getInfo(self):
        """ Returns info about junction."""
        out = 'b = {0:n}, Q_1 = {1:n}, rho = {2:n}'.format(self.b, self.getQ1(), self.getRho())
        return out

    def getQ1(self):
        """ Returns Q_1, the quality factor of the frequency dependent branch of the junction."""
        return math.sqrt(self.b)/(self.d + 1)

    def getRho(self):
        """ Returns (R*C)/(R_s*C_s) the quotient of the different junction time constants."""
        return self.e*self.b

    def getType(self):
        """ Returns type of junction."""
        return 'freq dep'

    def dx(self, x, t):
        """ Returns derivative for diff eq."""
        i_ = self.i
        b_ = self.b
        d_ = self.d
        e_ = self.e
        dp = x[1]
        dv = (i_ + i_n - x[1] - math.sin(x[0]))/b_
        dvc = e_*(v - vc)
        return [dp, dv, dvc]
    
    def applyI(self, i = 0.0, T=1000):
        """ Applies bias current to the junction for duration T.

            i: applied bias current
            T: duration of counting

            sets the phase to the most recent phase
            sets voltage to average over last three fifths 
            of T returns the average voltage."""
        self.i = i
        iT = int(T/self.dt) # renormalized time to integer values
        p = self.phase
        v = self.volt
        vc = self.v_c
        sumv = 0.0
        vtot = 0
        for t in xrange(iT):
            pvv = euler([p, v, vc], self.dx, self.dt)
            p, v, vc = pvv
            if t > 3*iT/5:
                sumv+=v
                vtot+=1
        self.setP(p)
        self.volt, self.v_c = v, vc
        return sumv/vtot

    def getPhaseVolt(self, i = 0.0):
        """ Applies current and returns phase and voltage at point."""
        self.i = i
        pvv = euler([self.phase, self.volt, self.v_c], self.dx, self.dt)
        self.phase, self.volt, self.v_c = pvv
        return pvv[:2]

class JJnFreq(JJFreq, JJn):
    """ A noisy Josephson Junction with frequency dependent circuit elements."""
    def __init__(self, b_c = 1.0, p = 0.0, v = 0.0, vc = 0.0, Q1 = 1.0, rho = 1.0, dt = .01, temp = 0.0):
        """ Initiates the junction.

            b_c: Stewart-McCumber damping constant
            p: the initial phase difference across the junction
            v: average voltage across the junction
            dt: timestep
            V_c: average voltage across the capacitor
            Q1: quality factor of the frequency dependent branch of the junction
            rho: quotient of the different junction time constants
            temp: unit-less temperature """

        JJFreq.__init__(self, b_c, p, v, vc, Q1, rho, dt)
        self.setTemp(temp)

    def getInfo(self, dt = .01):
        """ Returns info about junction."""
        out = 'b = {0:n}, temp = {1:n}, Q_1 = {2:n}, rho = {3:n}'.format(self.b, self.getTemp(), self.getQ1(), self.getRho())
        return out

    def getType(self):
        """ Returns type of junction."""
        return 'noisy freq dep'

    def getTemp(self):
        """ Returns normalized temperature at which the junction is operating."""
        return self.sig*self.sig*.5*self.dt

    def dx(self, x, t):
        """ Returns derivative for diff eq."""
        i_ = self.i
        b_ = self.b
        d_ = self.d
        e_ = self.e
        i_n = random.gauss(0, self.sig)
        i_n1 = random.gauss(0, self.sig * math.sqrt(d_))
        dp = x[1]
        dv = (i_ + i_n + i_n1 - x[1] - math.sin(x[0]) - d_*(x[1]-x[2]))/b_
        dvc = e_*(x[1] - x[2] - i_n1/d_)
        return [dp, dv, dvc]

class JJArray():
    """ A Josephson Junction Array"""
    def __init__(self, rows = 1, columns = 1, B_ = [[1.0]], iex = [0.0], dt = .01):
        """ Initiates the junction array.

            rows: number of rows of grains in the array
            columns : number of columns of grains in the array
            B_: damping matrix. It must be a square matrix of floats
            iex: external current multiplier. Multiplies external current on the grain by this (0.0 means no external current). All the positive components should add to 1 and the negatives to -1
            dt: timestep """
        self.cols = columns
        self.nodes = columns*rows
        self.pv = [0.0]*(self.nodes*2)
        self.B = B_
        self.Bi = matInv(B_)
        self.Iex = iex
        self.dt = float(dt)
        self.i = 0

    def getInfo(self):
        """ Returns info about junction array."""
        out = 'B = {0}, Iex = {1}'.format(self.B, self.Iex)
        return out

    def getType(self):
        """ Returns type of junction array."""
        return 'array'

    def applyI(self, i = 0.0, T=1000):
        """ Applies bias current to the junction array for duration T.

            i: applied bias current
            dt: timestep for numerics
            T: duration of counting

            sets the phase to the most recent phase 
            T returns the average voltage between the first 
            and last node."""
        self.i = i
        iT = int(T/self.dt) # renormalized time to integer values
        pv = self.pv
        sumv = 0.0
        vtot = 0
        for t in xrange(iT):
            # print(t)
            # print(pv[0] - pv[1])
            # print(pv[2] - pv[3])
            pv = euler(pv, self.dx, self.dt)
            if t > 3*iT/5:
                sumv += pv[self.nodes] - pv[self.nodes*2 - 1] # get voltage between first node and last node
                # print('{0},   {1}'.format(pv[self.nodes], pv[-1]))
                vtot += 1
        # print(pv)
        for i in xrange(self.nodes):
            pv[i] = math.fmod(pv[i], 2*math.pi)
        self.pv = pv
        return sumv/vtot

    def dx(self, x, t):
        """ Returns derivative for diff eq."""
        i_ = self.i
        N = self.nodes
        c = self.cols
        # print('c = {0}'.format(c))

        sm = [[0 for i in xrange(N)] for j in xrange(N)]
        # print('sm = {0}'.format(sm))
        for g in xrange(N):
            # print(sm[0][0])
            if g+c<N:
                sm[g][g+c] = x[N+g]-x[N+g+c] + math.sin(x[g]-x[g+c]) #+ random.gauss(0, sig)
                # print('1st if: {0}'.format(sm[0][0]))
                sm[g+c][g] = -1*sm[g][g+c]
                # print('1st if: {0}'.format(sm[0][0]))
            if (g+1)%c != 0:
                sm[g][g+1] = x[N+g]-x[N+g+1] + math.sin(x[g]-x[g+1]) #+ random.gauss(0, sig)
                sm[g+1][g] = -1*sm[g][g+1]
            # print(sm[0][0])

        # print('sm = {0}'.format(sm))
        s = [0]*N
        for k in xrange(N):
            s[k] = self.Iex[k]*i_ - sum(sm[k]) 

        # print('s = {0}'.format(s))
        dv = matmult(self.Bi, s)

        # print('dv = {0}'.format(dv))
        dx_ = [0.0]*(N*2)
        for k in xrange(N):
            dx_[k] = x[N+k]
            dx_[k+N] = dv[k]

        return dx_

    def getTemp(self):
        """ Returns normalized temperature at which the junction array is operating."""
        return 0