# -*- coding: utf-8 -*-

# This code is part of Schroedinger Newton Inhibitions of Dispersion
# (C) Copyright Andre Grossardt 2010-2023
# https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
# 
# This code is licensed under the MIT License (see LICENSE.txt for details)

"""
This script writes a file with tabulator separated data columns
t1 r90(t1)
t2 r90(t2)
...
where r90 is the radius within which 90% of the probability density is
contained.

Called with the following arguments:

r90.py path n max_t save_every w m dr dt file

    path: path to the data files
    n: number of grid points
    max_t: maximum time up to which the data is read
    save_every: save every nth time step (parameter of simulation)
    w: width of the initial gaussian
    m: mass of the particle
    dr: grid step size
    dt: time step size
    file: output file name
"""

import sys
import os
import numpy as np
from helpers import read_wave


def rc(psi, r, dr):
    """
    Calculate the radius within which 90% of the probability density is
    contained.

    Args:
        psi: wave function
        r: radial grid (array)
        dr: grid step size

    Returns:
        radius within which 90% of the probability density is contained
    """
    rho = (np.abs(psi) * r)**2
    rho /= rho.sum()
    i = 0
    s = rho[i]
    while ( s < .9 ):
        i += 1
        if i < len(rho):
            s += rho[i]
        else:
            s = 1.0
    return i * dr


if __name__ == "__main__":
    # arguments: path, n, max_t, save_every, w, m, dr, dt, file
    if len(sys.argv) < 9:
        print('need args: path, n, max_t, save_every, w, m, dr, dt')
        exit()
    path = sys.argv[1]
    if path[-1] != '/':
        path = path + '/'
    outpath = path
    path += 'data/'
    outpath += 'r90.dat'
    if not os.path.exists(os.path.dirname(path)):
        print('Error: Path does not exist')
        exit()
    n = int(sys.argv[2])
    max_t = int(sys.argv[3])
    save_every = int(sys.argv[4])
    w = float(sys.argv[5])
    m = float(sys.argv[6])
    dr = float(sys.argv[7])
    dt = float(sys.argv[8])
    
    if os.path.exists(outpath):
        print('Error: Outfile exists')
        exit()
    outfile = open(outpath,'w')
    outfile.write( "%e\t%e\n" % (0., w * 1.76796332416) )

    x = np.arange(n) * dr
    saved = max_t / save_every
    r = range(saved)
    t = (np.array(r) + 1) * save_every
    if max_t % save_every > 0:
        r.append(saved)
        t.append(max_t)
        saved += 1
    for i in r:
        time = (i+1) * save_every * dt * 1e-9
        print("t = %e s" % time)
        psi = read_wave(path, t[i], save_every)
        outfile.write( "%e\t%e\n" % (time, rc(psi, x, dr)) )
    outfile.close()