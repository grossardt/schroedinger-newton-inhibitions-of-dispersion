# -*- coding: utf-8 -*-

# This code is part of Schroedinger Newton Inhibitions of Dispersion
# (C) Copyright Andre Grossardt 2010-2023
# https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
# 
# This code is licensed under the MIT License (see LICENSE.txt for details)

"""
This script writes a select number of wave function(s) from the data files
to text files for further processing. The text files contain the radial
coordinate r and probability density rho(r) = |psi(r)|^2 as tabulator
separated columns:
r1 rho(r1)
r2 rho(r2)
...

Called with the following arguments:

plotfile.py path outfile-prefix n save_every dr timesteps

    path: path to the data files
    outfile-prefix: prefix of the output files
    n: number of grid points
    save_every: save every nth time step (parameter of simulation)
    dr: grid step size
    timesteps: time steps to plot (separated by spaces)

For every time step a separate file is created with the time value appended
to the file name as specified by outfile-prefix.
"""

import sys
import numpy as np
from helpers import read_wave, plot_colors



if __name__ == "__main__":
    # arguments: path, n, save_every, times
    if len(sys.argv) < 5:
        print('need arguments: ', end='')
        print('path, outfile-prefix, n, save_every, dr, timesteps')
    path = sys.argv[1]
    if path[-1] == '/':
        path = path[:-1]
    outfile_prefix = sys.argv[2]
    n = int(sys.argv[3])
    save_every = int(sys.argv[4])
    dr = float(sys.argv[5])
    
    times = []
    for arg in sys.argv[6:]:
        times.append(int(arg))
    
    if len(times) < 2:
        colors = ['blue']
    elif len(times) < 3:
        colors = ['red','blue']
    else:
        colors = plot_colors()
    
    r = np.arange(n) * dr;
    
    for i in range(len(times)):
        t = times[i]
        psi = read_wave(path, t, save_every)
        rho = (np.abs(psi) * r)**2
        outfile = open( ( "%s%d" % (outfile_prefix, t) ),'w')
        for j in range(len(rho)):
            outfile.write( "%e\t%e\n" % (r[j], rho[j]) )
        outfile.close()