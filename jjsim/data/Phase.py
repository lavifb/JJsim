import math
import JJs
import FileSetup as FS
import time, datetime

def phasePorts(js, i=.01, T=1000, mod = 100000, fl='test.dat'):
    """ Returns the phase portraits for the junctions.

        js: array of junctions
        i: bias current
        dt: timestep
        T: duration of junction averaging
        fl: data file to write """ 
    
    start_time = time.time()
    nw = datetime.datetime.now()
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    # write heading
    f.write('Phase Portrait \n')
    f.write('No. of Junction:    {0} \n'.format(len(js)))
    f.write('Type of Junctions:  {0} \n'.format(js[0].getType()))
    f.write('Junctions Info:     {0} \n'.format(js[0].getInfo()))
    f.write('dt, T:              {0}, {1} \n'.format(js[0].dt, T))
    f.write('i:                  {0} \n'.format(i))
    f.write('Date:               {0}/{1}/{2} \n'.format(nw.month, nw.day, nw.year))
    f.write('Time:               {0}:{1} \n'.format(nw.hour, nw.minute))

    out = 'Time    Phase    Voltage \n(t)    (p)    (v) \n'
    print('starting')
    for j in js:
        t = 0
        while t*j.dt<T:
            pv = j.getPhaseVolt(i)
            if t%mod == 0:
                out += '{0:.3f}    {1:.8f}    {2:.8f} \n'.format(t*j.dt, pv[0], pv[1])
            if pv[0] > 50:
                break
            t+=1
        out += '\n'
    f.write('Runtime:            {0} \n\n'.format(time.time() - start_time))
    f.write(out)
    f.close()
    print('done {0}'.format(fl))

def phasePort(j, i=.01, T=1000, mod = 100000, fl='test.dat'):
    """ Returns the phase portraits for a junction.

        j: junctions
        i: bias current
        dt: timestep
        T: duration of junction averaging
        fl: data file to write """ 
        start_time = time.time()
    nw = datetime.datetime.now()
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    # write heading
    f.write('Phase Portrait \n')
    f.write('No. of Junction:    {0} \n'.format(len(js)))
    f.write('Type of Junctions:  {0} \n'.format(js[0].getType()))
    f.write('Junctions Info:     {0} \n'.format(js[0].getInfo()))
    f.write('dt, T:              {0}, {1} \n'.format(js[0].dt, T))
    f.write('i:                  {0} \n'.format(i))
    f.write('Date:               {0}/{1}/{2} \n'.format(nw.month, nw.day, nw.year))
    f.write('Time:               {0}:{1} \n'.format(nw.hour, nw.minute))

    out = 'Time    Phase    Voltage \n(t)    (p)    (v) \n'
    print('starting')
    t = 0
    while t*j.dt<T:
        pv = j.getPhaseVolt(i)
        if t%mod == 0:
            out += '{0:.3f}    {1:.8f}    {2:.8f} \n'.format(t*j.dt, pv[0], pv[1])
        if pv[0] > 50:
            break
        t+=1
    f.write('Runtime:            {0} \n\n'.format(time.time() - start_time))
    f.write(out)
    f.close()
    print('done {0}'.format(fl))