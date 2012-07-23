import math
from jjsim import JJs
import FileSetup as FS
import time, datetime

def IVPlot(js, T, di = .01, i0 = 0, imax=1.5, fl='test.dat'):
    """ Produces data file with an IV curve that is averaged over the junctions.

        js: list of junctions used
        dt: time step 
        T: duration of junction averaging
        di: current step
        i0: initial current
        imax: max current
        fl: data file to write"""
    start_time = time.time()
    nw = datetime.datetime.now()
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    # write heading
    f.write('IV Plot \n')
    f.write('No. of Junctions:    {0} \n'.format(len(js)))
    f.write('Type of Junctions:  {0} \n'.format(js[0].getType()))
    f.write('Junctions Info:     {0} \n'.format(js[0].getInfo()))
    f.write('dt, T:              {0}, {1} \n'.format(js[0].dt, T))
    f.write('current range, di:  {1}-{2}, {0} \n'.format(di, i0, imax))
    f.write('Date:               {0}/{1}/{2} \n'.format(nw.month, nw.day, nw.year))
    f.write('Time:               {0}:{1} \n'.format(nw.hour, nw.minute))

    i = i0
    out = 'Current    Voltage \n(i)        (v) \n'
    while (i < imax):
        out = ''.join([out, '{0:.5f}'.format(i)])
        sumv = 0.0
        vtot = 0
        for j in js:
            v = j.applyI(i, T)
            sumv += v
            vtot += 1
        out = ''.join([out, '    {0:+.8f} \n'.format(sumv/vtot)])
        print('{0} run for {1} seconds'.format(i, time.time() - start_time))
        i += di
    f.write('Runtime:            {0} \n\n'.format(time.time() - start_time))
    f.write(out)
    f.close()
    print('done with {0}'.format(fl))

def hyst(js, T, di = .01, i0 = 0, imax=1.5, fl='test.dat'):
    """ Produces data file with an hysteric IV curve that is averaged over the junctions.

        js: list of junctions used
        dt: time step 
        T: duration of junction averaging
        di: current step
        i0: initial current
        imax: max current
        fl: data file to write"""
    start_time = time.time()
    nw = datetime.datetime.now()
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    # write heading
    f.write('Hysteric IV Plot \n')
    f.write('No. of Junctions:    {0} \n'.format(len(js)))
    f.write('Type of Junctions:  {0} \n'.format(js[0].getType()))
    f.write('Junctions Info:     {0} \n'.format(js[0].getInfo()))
    f.write('dt, T:              {0}, {1} \n'.format(js[0].dt, T))
    f.write('current range, di:  {1}-{2}, {0} \n'.format(di, i0, imax))
    f.write('Date:               {0}/{1}/{2} \n'.format(nw.month, nw.day, nw.year))
    f.write('Time:               {0}:{1} \n'.format(nw.hour, nw.minute))

    i = i0
    out = 'Current    Voltage \n(i)        (v) \n'
    while (i < imax):
        out = ''.join([out, '{0:.5f}'.format(i)])
        sumv = 0.0
        vtot = 0
        for j in js:
            v = j.applyI(i, T)
            sumv += v
            vtot += 1
        out = ''.join([out, '    {0:+.8f} \n'.format(sumv/vtot)])
        print('{0} run for {1} seconds'.format(i, time.time() - start_time))
        i += di
    i = imax
    print('Lowering current')
    while (i > i0):
        out = ''.join([out, '{0:.5f}'.format(i)])
        sumv = 0.0
        vtot = 0
        for j in js:
            v = j.applyI(i, T)
            sumv += v
            vtot += 1
        out = ''.join([out, '    {0:+.8f} \n'.format(sumv/vtot)])
        print('{0} run for {1} seconds'.format(i, time.time() - start_time))
        i -= di
    f.write('Runtime:            {0} \n\n'.format(time.time() - start_time))
    f.write(out)
    f.close()
    print('done  {0}'.format(fl))

def allIVPlot(js, T = 1000, di = .01, i0 = 0.0, imax=1.5, fl='test.dat'):
    """ Produces data file with an IV curve that is averaged over the junctions
        along with data from each individual junction.

        js: list of junctions used
        dt: time step 
        T: duration of junction averaging
        di: current step
        i0: initial current
        imax: max current
        fl: data file to write"""
    start_time = time.time()
    nw = datetime.datetime.now()
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    # write heading
    f.write('IV Plot with individual junctions\n')
    f.write('No. of Junctions:    {0} \n'.format(len(js)))
    f.write('Type of Junctions:  {0} \n'.format(js[0].getType()))
    f.write('Junctions Info:     {0} \n'.format(js[0].getInfo()))
    f.write('dt, T:              {0}, {1} \n'.format(js[0].dt, T))
    f.write('current range, di:  {1}-{2}, {0} \n'.format(di, i0, imax))
    f.write('Date:               {0}/{1}/{2} \n'.format(nw.month, nw.day, nw.year))
    f.write('Time:               {0}:{1} \n'.format(nw.hour, nw.minute))

    i = i0
    out = 'Current    Voltage \n(i)        (v) \n'
    while (i < imax):
        out = ''.join([out, '{0:.5f}'.format(i)])
        vs = []
        sumv = 0.0
        vtot = 0
        for j in js:
            v = j.applyI(i, T)
            vs.append(v)
            sumv += v
            vtot += 1
        out = ''.join([out, '    {0:+.8f} \n'.format(sumv/vtot)])
        for v_ in vs:
            out = ''.join([out, '    {0:+.8f}'.format(v_)])
        out = ''.join([out, '\n'])
        print('{0} run for {1} seconds'.format(i, time.time() - start_time))
        i += di
    f.write('Runtime:            {0} \n\n'.format(time.time() - start_time))
    f.write(out)
    f.close()
    print('done with {0}'.format(fl))

def allHyst(js, T, di = .01, i0 = 0, imax=1.5, fl='test.dat'):
    """ Produces data file with an hysteric IV curve that is averaged over the junctions
        along with data from each individual junction.

        js: list of junctions used
        dt: time step 
        T: duration of junction averaging
        di: current step
        i0: initial current
        imax: max current
        fl: data file to write"""
    start_time = time.time()
    nw = datetime.datetime.now()
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    # write heading
    f.write('Hysteric IV Plot with individual junctions\n')
    f.write('No. of Junctions:    {0} \n'.format(len(js)))
    f.write('Type of Junctions:  {0} \n'.format(js[0].getType()))
    f.write('Junctions Info:     {0} \n'.format(js[0].getInfo()))
    f.write('dt, T:              {0}, {1} \n'.format(js[0].dt, T))
    f.write('current range, di:  {1}-{2}, {0} \n'.format(di, i0, imax))
    f.write('Date:               {0}/{1}/{2} \n'.format(nw.month, nw.day, nw.year))
    f.write('Time:               {0}:{1} \n'.format(nw.hour, nw.minute))

    i = i0
    out = 'Current    Voltage \n(i)        (v) \n'
    while (i < imax):
        out = ''.join([out, '{0:.5f}'.format(i)])
        vs = []
        sumv = 0.0
        vtot = 0
        for j in js:
            v = j.applyI(i, T)
            vs.append(v)
            sumv += v
            vtot += 1
        out = ''.join([out, '    {0:+.8f} \n'.format(sumv/vtot)])
        for v_ in vs:
            out = ''.join([out, '    {0:+.8f}'.format(v_)])
        out = ''.join([out, '\n'])
        print('{0} run for {1} seconds'.format(i, time.time() - start_time))
        i += di
    i = imax
    print('Lowering current')
    while (i > i0):
        out = ''.join([out, '{0:.5f}'.format(i)])
        vs = []
        sumv = 0.0
        vtot = 0
        for j in js:
            v = j.applyI(i, T)
            vs.append(v)
            sumv += v
            vtot += 1
        out = ''.join([out, '    {0:+.8f} \n'.format(sumv/vtot)])
        for v_ in vs:
            out = ''.join([out, '    {0:+.8f}'.format(v_)])
        out = ''.join([out, '\n'])
        print('{0} run for {1} seconds'.format(i, time.time() - start_time))
        i -= di
    f.write('Runtime:            {0} \n\n'.format(time.time() - start_time))
    f.write(out)
    f.close()
    print('done  {0}'.format(fl))
