import os

def fileSetup(fl='test.dat'):
    """ Sets up the file. Returns file to be written."""
    dir = os.path.dirname(fl)
    print(dir)
    if not os.path.exists(dir) and dir != "":
        os.makedirs(dir)
    fil = fl
    k = 2
    while os.path.exists(fil):
        fil = '{0}_{1}.dat'.format(fl[:-4], k)
        k+=1
    print(fil)
    return fil