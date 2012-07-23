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
    
def dot(v1, v2):
    'dot product of two vectors.'
    sum = 0
    for i in xrange(len(v1)):
        sum += v1[i] * v2[i]
    return sum

def matmult(m, v):
    'multiply matrix by a vector.'
    return [dot(r, v) for r in m]

def gauss_jordan(m, eps = 1.0/(10**10)):
    """Puts given matrix (2D array) into the Reduced Row Echelon Form.
        Returns True if successful, False if 'm' is singular."""
    (h, w) = (len(m), len(m[0]))
    for y in range(0,h):
        maxrow = y
        for y2 in range(y+1, h):    # Find max pivot
            if abs(m[y2][y]) > abs(m[maxrow][y]):
              maxrow = y2
        (m[y], m[maxrow]) = (m[maxrow], m[y])
        if abs(m[y][y]) <= eps:     # Singular?
            return False
        for y2 in range(y+1, h):    # Eliminate column y
            c = m[y2][y] / m[y][y]
            for x in range(y, w):
                m[y2][x] -= m[y][x] * c
    for y in range(h-1, 0-1, -1): # Backsubstitute
        c  = m[y][y]
        for y2 in range(0,y):
            for x in range(w-1, y-1, -1):
                m[y2][x] -=  m[y][x] * m[y2][y] / c
        m[y][y] /= c
        for x in range(h, w):       # Normalize row y
            m[y][x] /= c
    return True

def matInv(M):
    """ Returns the inverse of the matrix M."""
    m2 = [row[:]+[int(i==j) for j in range(len(M) )] for i,row in enumerate(M) ]
    return [row[len(M[0]):] for row in m2] if gauss_jordan(m2) else None