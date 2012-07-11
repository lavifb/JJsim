import math
import JJs
import FileSetup as FS
import time, datetime

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