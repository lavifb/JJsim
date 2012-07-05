import math
import RungeKutta as RK
import random

class JJ:
    """ A Josephson Junction"""
    def __init__(self, b_c = 1.0, p = 0.0, v = 0.0, dt = .01):
        """ Initiates the junction.

            b_c : Stewart-McCumber damping constant
            p: the initial phase difference across the junction
            v: average voltage across the junction
            dt: timestep """
        self.phase = float(p)
        self.volt = float(v)
        self.b = float(b_c)
        # b>1 is underdmaped. b<1 is overdamped
        self.dt = float(dt)

    def getInfo(self):
        """ Returns info about junction."""
        out = 'b = {0}'.format(self.b)
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
            sets voltage to average over last fifth of T
            returns the average voltage."""
        v = self.volt
        t = 0
        sumv = 0.0
        vtot = 0
        while (t*self.dt<T):
            pv = RK.rk4_2(self.phase, v, self.getA(i), self.dt)
            self.setP(pv[0])
            v = pv[1]
            if t > 4*T/5:
                sumv+=v
                vtot+=1
            t += 1
        self.volt = sumv/vtot
        return self.volt

    def getA(self, i):
        """ Returns the acceleration function to use in RK4."""
        b_ = self.b
        def a(x, v, dt):
            return (i - v - math.sin(x))/b_
        return a

    def setP(self, p):
        """ Sets the phase, taking it mod 2*pi"""
        self.phase = math.fmod(p, 2*math.pi)
        pass

    def getTemp(self):
        """ Returns normalized temperature at which the junction is operating."""
        return 0

    def getPhaseVolt(self, i = 0.0):
        """ Applies current and returns phase and voltage at point."""
        v = self.volt
        t = 0
        pv = RK.rk4_2(self.phase, v, self.getA(i), self.dt)
        self.phase = pv[0]
        self.volt = pv[1]
        return pv

class JJn(JJ):
    """A Josephson Junction with thermal noise"""
    def __init__(self, b_c = 1.0, p = 0.0, v = 0.0, dt = .01, temp = 0.0):
        """ Initiates the junction.

            b_c : Stewart-McCumber damping constant
            p: the initial phase difference across the junction
            v: average voltage across the junction
            dt: timestep 
            temp: unit-less temperature """

        JJ.__init__(self, b_c, p, v, dt)
        self.setTemp(temp)

    def getInfo(self):
        """ Returns info about junction."""
        out = 'b = {0}, temp = {1}'.format(self.b, self.getTemp())
        return out

    def getType(self):
        """ Returns type of junction"""
        return 'noisy'

    def getA(self, i):
        """ Returns the acceleration function to use in RK4."""
        b_ = self.b
        i_n = random.gauss(0, self.sig)
        def a(x, v, dt):
            return (i + i_n - v - math.sin(x))/b_
        return a

    def setSig(self, s):
        """ Sets the st. dev. for thermal noise current."""
        self.sig = s
        pass

    def getTemp(self):
        """ Returns normalized temperature at which the junction is operating."""
        return self.sig*self.sig*.5*self.dt

    def setTemp(self, temp):
        """ Sets self.sig to the appropriate value for the temperature."""
        self.sig = math.sqrt((temp + temp)/self.dt)
        pass

class JJFreq(JJ):
    """ A Josephson Junction with frequency dependent circuit elements."""
    def __init__(self, b_c = 1.0, p = 0.0, v = 0.0, dt = .01, V_c = 0.0, D = 1.0, E = 1.0):
        """ Initiates the junction.

            b_c : Stewart-McCumber damping constant
            p: the initial phase difference across the junction
            v: average voltage across the junction
            dt: timestep 
            V_c: average voltage across the capacitor
            D: (Q_0/Q_1 - 1) of the configuration
            E: p/(Q_0^2) = 1/(R*C*R_s*C_s*(w_p)^2) of the configuration """

        JJ.__init__(self, b_c, p, v, dt)
        self.v_c = float(V_c)   # voltage across shunting capacitor
        self.d = float(D)       # (Q_0/Q_1 -1)
        self.e = float(E)       # p/(Q_0^2) = 1/(R*C*R_s*C_s*(w_p)^2)

    def getInfo(self):
        """ Returns info about junction."""
        out = 'b = {0}, d = {1}, e = {2}'.format(self.b, self.d, self.e)
        return out

    def getType(self):
        """ Returns type of junction."""
        return 'freq dep'

    def getDp(self):
        """ Returns derivative of phase to use in RK4."""
        def dp(p, v, vc, t):
            return v
        return dp

    def getDv(self, i):
        """ Returns derivative of voltage to use in RK4."""
        b_ = self.b
        d_ = self.d
        def dv(p, v, vc, t):
            return (i - v - math.sin(p) - d_*(v-vc))/b_
        return dv

    def getDvc(self):
        """ Returns derivative of v_c to use in RK4."""
        e_ = self.e
        def dvc(p, v, vc, t):
            return e_*(v - vc)
        return dvc
    
    def applyI(self, i = 0.0, T=1000):
        """ Applies bias current to the junction for duration T.

            i: applied bias current
            T: duration of counting

            sets the phase to the most recent phase
            sets voltage to average over last fifth of T
            returns the average voltage."""
        v = self.volt
        vc = self.v_c
        t = 0
        sumv = 0.0
        vtot = 0
        while (t*self.dt<T):
            pvv = RK.rk4_3(self.phase, v, vc, self.getDp(), self.getDv(i), self.getDvc(), self.dt)
            self.setP(pvv[0])
            v = pvv[1]
            vc = pvv[2]
            if t > 4*T/5:
                sumv+=v
                vtot+=1
            t += 1
        self.volt = sumv/vtot
        self.v_c = vc
        return self.volt

class JJnFreq(JJFreq, JJn):
    """ A noisy Josephson Junction with frequency dependent circuit elements."""
    def __init__(self, b_c = 1.0, p = 0.0, v = 0.0, dt = .01, V_c = 0.0, D = 1.0, E = 1.0, temp = 0.0):
        """ Initiates the junction.

            b_c : Stewart-McCumber damping constant
            p: the initial phase difference across the junction
            v: average voltage across the junction
            dt: timestep
            V_c: average voltage across the capacitor
            D: (Q_0/Q_1 - 1) of the configuration
            E: p/(Q_0^2) = 1/(R*C*R_s*C_s*(w_p)^2) of the configuration
            temp: unit-less temperature """

        JJFreq.__init__(self, b_c, p, v, dt, V_c, D, E)
        self.setTemp(temp, dt)

    def getInfo(self, dt = .01):
        """ Returns info about junction."""
        out = 'b = {0}, temp = {1}, d = {2}, e = {3}'.format(self.b, self.getTemp(dt), self.d, self.e)
        return out

    def getType(self):
        """ Returns type of junction."""
        return 'noisy freq dep'

    def getDv(self, i):
        """ Returns derivative of voltage to use in RK4."""
        b_ = self.b
        d_ = self.d
        i_n = random.gauss(0, self.sig)
        i_n1 = random.gauss(0, self.sig * math.sqrt(d_))
        def dv(p, v, vc, t):
            return (i + i_n + i_n1 - v - math.sin(p) - d_*(v-vc))/b_
        return dv

    def getDvc(self):
        """ Returns derivative of v_c to use in RK4."""
        d_ = self.d
        e_ = self.e
        i_n1 = random.gauss(0, self.sig * math.sqrt(d_))
        def dvc(p, v, vc, t):
            return e_*(v - vc - i_n1/d_)
        return dvc
