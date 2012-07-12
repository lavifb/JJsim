def RK4(xs, fxs, dt, t = 0):
    """ Returns final xs tuple after time dt has passed.

    xs: initial conditions (tuple)
    fxs: dx/dt functions, fx(xs, t) (tuple: elements must be callable)
    dt: timestep (number)
    t: initial time (number)"""

    l = len(xs)
    if l != len(fxs):
        raise Exception('Need as many functions as variables')

    fx1 = []
    for fx in fxs:
        fx1.append(dt * fx(xs, t))

    x2 = []
    for i in xrange(0,l):
        x2.append(xs[i] + fx1[i]*0.5)
    fx2 = []
    for fx in fxs:
        fx2.append(dt * fx(x2, t + dt*0.5))

    x3 = []
    for i in xrange(0,l):
        x3.append(xs[i] + fx2[i]*0.5)
    fx3 = []
    for fx in fxs:
        fx3.append(dt * fx(x3, t + dt*0.5))

    x4 = []
    for i in xrange(0,l):
        x4.append(xs[i] + fx3[i])
    fx4 = []
    for fx in fxs:
        fx3.append(dt * fx(x4, t + dt))

    xf = []
    for i in xrange(0,l):
        xf.append(xs[i] + (fx1[i] + fx2[i] + fx2[i] + fx3[i] + fx3[i] + fx4[i])/6.0)

    return xf

def rk4_1(t, x, f, dt):
    """ Returns final position after time dt has passed.

    t: time (number)
    x: position (number)
    f: derivative function f(x,t) (must be callable)
    dt: timestep (number)"""

    k1 = dt * f(x, t)
    k2 = dt * f(x + 0.5*k1, t + 0.5*dt)
    k3 = dt * f(x + 0.5*k2, t + 0.5*dt)
    k4 = dt * f(x + k3, t + dt)

    xf = x + (k1 + k2 + k2 + k3 + k3 + k4)/6.0
    
    return xf

def euler(x0, f, dt, t = 0):
    xf = x0
    dx = f(x0, t)
    for i in xrange(len(x0)):
        xf[i] = x0[i] + dt*dx[i]
    return xf