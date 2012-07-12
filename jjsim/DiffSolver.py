def rk4(p, a, dt):
    """Returns final (position, velocity) tuple after time dt has passed.

    p: point in form numpy.array([x,v])
    a: acceleration function a(x,v,dt) (must be callable)
    dt: timestep (number)"""
    x1 = p[0]
    v1 = p[1]
    a1 = a(x1, v1, 0)

    x2 = p[0] + 0.5*v1*dt
    v2 = p[1] + 0.5*a1*dt
    a2 = a(x2, v2, dt/2.0)

    x3 = p[0] + 0.5*v2*dt
    v3 = p[1] + 0.5*a2*dt
    a3 = a(x3, v3, dt/2.0)

    x4 = p[0] + v3*dt
    v4 = p[1] + a3*dt
    a4 = a(x4, v4, dt)

    xf = p[0] + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
    vf = p[1] + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)

    pf = [xf, vf]
    return pf

def rk4_2(x, v, a, dt):
    """ Returns final (position, velocity) tuple after time dt has passed.

    x: position (number)
    v: velocity (number)
    a: acceleration function a(x,v,dt) (must be callable)
    dt: timestep (number)"""
    x1 = x
    v1 = v
    a1 = a(x1, v1, 0)

    x2 = x + 0.5*v1*dt
    v2 = v + 0.5*a1*dt
    a2 = a(x2, v2, dt/2.0)

    x3 = x + 0.5*v2*dt
    v3 = v + 0.5*a2*dt
    a3 = a(x3, v3, dt/2.0)

    x4 = x + v3*dt
    v4 = v + a3*dt
    a4 = a(x4, v4, dt)

    xf = x + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
    vf = v + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)

    pf = [xf, vf]
    return pf

def rk4_3(x, y, z, fx, fy, fz, dt, t = 0):
    """ Returns final (x, y, z) tuple after time dt has passed.

    x: first variable initial (number)
    y: second variable initial (number)
    z: third variable initial (number)
    fx: dx/dt function, fx(x,y,z,t) (must be callable)
    fy: dy/dt function, fy(x,y,z,t) (must be callable)
    fz: dz/dt function, fz(x,y,z,t) (must be callable)
    dt: timestep (number)
    t: initial time (number)"""

    fx1 = dt * fx(x, y, z, t)
    fy1 = dt * fy(x, y, z, t)
    fz1 = dt * fz(x, y, z, t)

    fx2 = dt * fx(x + fx1*0.5, y + fy1*0.5, z + fz1*0.5, t + dt*0.5)
    fy2 = dt * fy(x + fx1*0.5, y + fy1*0.5, z + fz1*0.5, t + dt*0.5)
    fz2 = dt * fz(x + fx1*0.5, y + fy1*0.5, z + fz1*0.5, t + dt*0.5)

    fx3 = dt * fx(x + fx2*0.5, y + fy2*0.5, z + fz2*0.5, t + dt*0.5)
    fy3 = dt * fy(x + fx2*0.5, y + fy2*0.5, z + fz2*0.5, t + dt*0.5)
    fz3 = dt * fz(x + fx2*0.5, y + fy2*0.5, z + fz2*0.5, t + dt*0.5)

    fx4 = dt * fx(x + fx3, y + fy3, z + fz3, t + dt)
    fy4 = dt * fy(x + fx3, y + fy3, z + fz3, t + dt)
    fz4 = dt * fz(x + fx3, y + fy3, z + fz3, t + dt)

    xf = x + (fx1 + fx2 + fx2 + fx3 + fx3 + fx4)/6.0
    yf = y + (fy1 + fy2 + fy2 + fy3 + fy3 + fy4)/6.0
    zf = z + (fz1 + fz2 + fz2 + fz3 + fz3 + fz4)/6.0

    xyz = [xf, yf, zf]
    return xyz

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
    for i in range(0,l):
        x2.append(xs[i] + fx1[i]*0.5)
    fx2 = []
    for fx in fxs:
        fx2.append(dt * fx(x2, t + dt*0.5))

    x3 = []
    for i in range(0,l):
        x3.append(xs[i] + fx2[i]*0.5)
    fx3 = []
    for fx in fxs:
        fx3.append(dt * fx(x3, t + dt*0.5))

    x4 = []
    for i in range(0,l):
        x4.append(xs[i] + fx3[i])
    fx4 = []
    for fx in fxs:
        fx3.append(dt * fx(x4, t + dt))

    xf = []
    for i in range(0,l):
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