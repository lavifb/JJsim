import math
import JJs
import FileSetup as FS
import time, datetime

def dampingPlot(j, w0, wmax, estep):
	start_time = time.time()
    nw = datetime.datetime.now()
    fl = FS.fileSetup(fl)
    f = open(fl, 'w')
    # write heading
    f.write('Frequency Dependent Damping Plot \n')
    f.write('Type of Junction:     {0} \n'.format(js[0].getType()))
    f.write('Junctions Info:       {0} \n'.format(js[0].getInfo()))
    f.write('frequency range, dw:  {1}-{2}, {0} \n'.format(dw, w0, wmax))
    f.write('Date:                 {0}/{1}/{2} \n'.format(nw.month, nw.day, nw.year))
    f.write('Time:                 {0}:{1} \n'.format(nw.hour, nw.minute))

    w = w0
    out = 'Frequency    Damping \n(w)    (b) \n'
    while (w < wmax):
        out += '{0:.3e}'.format(w)
        b = b0((1 + math.sqrt(b0/b1)*b0*p**-2*(w/wp)**2)/(1 + b0**2/b1*p**-2*(w/wp)**2))**2
        out += '    {0:+.8e} \n'.format(b)
        w *= 10**estep
    f.write('Runtime:            {0} \n\n'.format(time.time() - start_time))
    f.write(out)
    f.close()
    print('done with {0}'.format(fl))