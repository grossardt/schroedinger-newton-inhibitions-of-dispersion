# -*- coding: utf-8 -*-

# This code is part of Schroedinger Newton Inhibitions of Dispersion
# (C) Copyright Andre Grossardt 2010-2023
# https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
# 
# This code is licensed under the MIT License (see LICENSE.txt for details)

"""
This script writes a file with tabulator separated data columns
t1 rmax(t1)
t2 rmax(t2)
...
where rmax is the peak of the radial probability density r^2 |psi|^2.

Called with the following arguments:

rmax.py path n max_t save_every w m dr dt file

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


def rmax(psi, r, dr):
    """
    Calculate the peak of the radial probability density r^2 |psi|^2
    for a given wave function psi.

    Args:
        psi: wave function
        r: radial grid (array)
        dr: grid step size

    Returns:
        peak of the radial probability density
    """
    rho = (np.abs(psi) * r)**2
    return rho.argmax() * dr


if __name__ == "__main__":
    # arguments: path, n, max_t, save_every, w, m, dr, dt, file
    if len(sys.argv) < 10:
        print('need args: path, n, max_t, save_every, w, m, dr, dt, file')
        exit()
    path = sys.argv[1]
    if path[-1] != '/':
        path = path + '/'
    path += 'data/'
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
    if os.path.exists(sys.argv[9]):
        print('Error: Outfile exists')
        exit()
    outfile = open(sys.argv[9],'w')
    outfile.write( "%e\t%e\n" % (0., w) )

    x = np.arange(n) * dr
    saved = max_t / save_every
    r = range(saved)
    t = (np.array(r) + 1) * save_every
    if max_t % save_every > 0:
        r.append(saved)
        t.append(max_t)
        saved += 1
    for i in r:
        time = (i+1) * save_every * dt
        psi = read_wave(path, n, save_every, t[i])
        outfile.write( "%e\t%e\n" % (time, rmax(psi, x, dr)) )
    outfile.close()