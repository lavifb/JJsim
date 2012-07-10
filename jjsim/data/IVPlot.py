import math
import JJs
import FileSetup as FS
import time, datetime

def IVPlotsDat(js, T, di = .01, i0 = 0, imax=1.5, fl='test.dat'):
    """ Produces data file with an IV curve that is averaged over the junctions
        along with the data from the actual junctions"""
    start_time = time.time()
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    f.write('# IV simulation data file\n')
    f.write('# This the average of {2} junctions with b = {0} and sig = {1}\n'.format(js[0].b, js[0].sig, len(js)))
    f.write('# The simulation is run with dt = {0} and T = {1}\n'.format(dt, T))
    f.write('# i       AveVolt         v1 ...\n')
    r = 0
    i = i0
    out = ''
    while (i < imax):
        out += '{0:.3f}'.format(i)
        vs = []
        sumv = 0.0
        vtot = 0
        for j in js:
            v = j.applyI(i, T)
            vs.append(v)
            sumv += v
            vtot += 1
        out += '  {0:+.8f}'.format(sumv/vtot)
        for v_ in vs:
            out += '  {0:+.8f}'.format(v_)
        out += '\n'
        print('{0} run for {1} seconds'.format(i, time.time() - start_time))
        i += di
    f.write('# Runtime: {0} seconds'.format(time.time() - start_time))
    f.write(out)
    f.close()
    print('done with {0}'.format(fl))

def hystDat(js, T, di = .01, i0 = 0, imax=1.5, fl='test.dat'):
    """ Produces data file with an hysteric IV curve that is averaged over the junctions
        along with the data from the actual junctions"""
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    f.write('# IV simulation data file\n')
    f.write('# This the average of {2} junctions with b = {0} and sig = {1}\n'.format(js[0].b, js[0].sig, len(js)))
    f.write('# The simulation is run with dt = {0} and T = {1}\n'.format(dt, T))
    f.write('# i       AveVolt         v1 ...\n')
    r = 0
    i = i0
    while (i < imax):
        f.write('{0:.3f}'.format(i))
        vs = []
        sumv = 0.0
        vtot = 0
        for j in js:
            v = j.applyI(i, T)
            vs.append(v)
            sumv += v
            vtot += 1
        f.write('  {0:+.8f}'.format(sumv/vtot))
        for v_ in vs:
            f.write('  {0:+.8f}'.format(v_))
        f.write('\n')
        i += di
    print('Going back down')
    i = imax
    while (i > i0):
        f.write('{0:.3f}'.format(i))
        vs = []
        sumv = 0.0
        vtot = 0
        for j in js:
            v = j.applyI(i, T)
            vs.append(v)
            sumv += v
            vtot += 1
        f.write('  {0:+.8f}'.format(sumv/vtot))
        for v_ in vs:
            f.write('  {0:+.8f}'.format(v_))
        f.write('\n')
        i-=di
    f.close()
    print('done  {0}'.format(fl))

def aveHystDat(js, T, di = .01, i0 = 0, imax=1.5, fl='test.dat'):
    """ Produces data file with an hysteric IV curve that is averaged over the junctions

        js: list of junctions used
        dt: time step 
        T: duration of junction averaging
        di: current step
        i0: initial current
        imax: max current
        fl: data file to write"""
    start_time = time.time()
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    f.write('# IV simulation data file\n')
    f.write('# This the average of {2} junctions with b = {0} and sig = {1}\n'.format(js[0].b, js[0].sig, len(js)))
    f.write('# The simulation is run with dt = {0} and T = {1}\n'.format(dt, T))
    f.write('# i       AveVolt\n')
    r = 0
    i = i0
    while (i < imax):
        f.write('{0:.3f}'.format(i))
        sumv = 0.0
        vtot = 0
        for j in js:
            sumv += j.applyI(i, T)
            vtot += 1
        f.write('  {0:+.8f}'.format(sumv/vtot))
        f.write('\n')
        print('{0} run for {1} seconds'.format(i, time.time() - start_time))
        i += di
    print('Going back down')
    i = imax
    while (i > i0):
        f.write('{0:.3f}'.format(i))
        sumv = 0.0
        vtot = 0
        for j in js:
            sumv += j.applyI(i, T)
            vtot += 1
        f.write('  {0:+.8f}'.format(sumv/vtot))
        f.write('\n')
        print('{0} run for {1} seconds'.format(i, time.time() - start_time))
        i-=di
    f.write('# Runtime: {0} seconds'.format(time.time() - start_time))
    f.close()
    print('done  {0}'.format(fl))

def aveIVPlot(js, T, di = .01, i0 = 0, imax=1.5, fl='test.dat'):
    """ Produces data file with an IV curve that is averaged over the junctions

        js: list of junctions used
        dt: time step 
        T: duration of junction averaging
        di: current step
        i0: initial current
        imax: max current
        fl: data file to write"""
    start_time = time.time()
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    f.write('# IV simulation data file\n')
    f.write('# This the average of {2} junctions with b = {0} and sig = {1}\n'.format(js[0].b, js[0].sig, len(js)))
    f.write('# The simulation is run with dt = {0} and T = {1}\n'.format(dt, T))
    f.write('# i       AveVolt\n')
    r = 0
    i = i0
    while (i < imax):
        f.write('{0:.3f}'.format(i))
        sumv = 0.0
        vtot = 0
        for j in js:
            sumv += j.applyI(i, T)
            vtot += 1
        f.write('  {0:+.8f}'.format(sumv/vtot))
        f.write('\n')
        print('{0} run for {1} seconds'.format(i, time.time() - start_time))
        i += di
    f.write('# Runtime: {0} seconds'.format(time.time() - start_time))
    f.close()
    print('done  {0}'.format(fl))

def phasePorts(js, i=.01, T=1000, fl='test.dat'):
    """ Returns the phase portraits for the junctions.

        js: array of junctions
        i: bias current
        dt: timestep
        T: duration of junction averaging
        fl: data file to write """ 
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    print('initializing junctions...')
    print('starting')
    count = 0
    for j in js:
        t = 0
        mod = 0
        while t<T:
            pv = j.getPhaseVolt(i)
            mod += 1
            if mod%100000 == 0:
                f.write('{0:.8f} {1:.8f} \n'.format(pv[0], pv[1]))
            if pv[0] > 50:
                break
            t+=j.dt
        count += 1
        f.write('\n')
    f.close()
    print('done {0}'.format(fl))

def phasePort(j, i=.01, T=1000, fl='test.dat'):
    """ Returns the phase portraits for a junction.

        j: junctions
        i: bias current
        dt: timestep
        T: duration of junction averaging
        fl: data file to write """ 
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    print('initializing junction...')
    print('starting')
    count = 0
    t = 0
    while t<T:
        pv = j.getPhaseVolt(i)
        f.write('{0:.8f} {1:.8f} \n'.format(pv[0], pv[1]))
        if pv[0] > 50:
            break
        t+=j.dt
    print('junction done'.format(count))
    count += 1
    f.close()
    print('done {0}'.format(fl))