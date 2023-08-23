# -*- coding: utf-8 -*-

# This code is part of Schroedinger Newton Inhibitions of Dispersion
# (C) Copyright Andre Grossardt 2010-2023
# https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
# 
# This code is licensed under the MIT License (see LICENSE.txt for details)

"""
This script prints the time for which the relative deviation of the
width of the wave function from the free solution is above a given
threshold.

Called with the following arguments:

masstime.py path n max_t save_every w m dr dt deviation

    path: path to the data files
    n: number of grid points
    max_t: maximum time up to which the data is read
    save_every: save every nth time step (parameter of simulation)
    w: width of the initial gaussian
    m: mass of the particle
    dr: grid step size
    dt: time step size
    deviation: desired relative deviation in width
"""

import sys
import os
import numpy as np
from helpers import free_solution, read_wave


def halfwidth(psi):
    """
    Calculate half the width of the wave function, i.e. the grid point
    index where the wave function goes below half the maximum value.

    Args:
        psi: wave function

    Returns:
        half the width of the wave function
    """
    rho = (np.abs(psi))**2
    m = rho.max()
    for i in range(len(rho)):
        if rho[i] < m/2.:
            return float(i)
    return float('inf')


if __name__ == "__main__":
    # arguments: path, n, max_t, save_every, w, m, dr, dt, rel. error
    if len(sys.argv) < 7:
        print('need args: path, n, max_t, save_every, w, m, dr, dt, rel. error')
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
    err = float(sys.argv[9])

    # determine time when width psi/free >= rel. error
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
        f = free_solution(w, m, time, n, dr)
        g = read_wave(path, t[i], save_every)
        if abs(1. - halfwidth(g)/halfwidth(f)) > err:
            print("%e \t %e" % (m, time / 1.0e9))
            exit()
    print("No difference to free solution")
