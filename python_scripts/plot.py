# -*- coding: utf-8 -*-

# This code is part of Schroedinger Newton Inhibitions of Dispersion
# (C) Copyright Andre Grossardt 2010-2023
# https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
# 
# This code is licensed under the MIT License (see LICENSE.txt for details)

"""
This script plots a (select number of) wave function(s) from the data files.

Called with the following arguments:

plot.py path n save_every times

    path: path to the data files
    n: number of grid points
    save_every: save every nth time step (parameter of simulation)
    times: time steps to plot (separated by spaces)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from helpers import read_wave, plot_colors


if __name__ == "__main__":
    # arguments: path, n, save_every, times
    if len(sys.argv) < 5:
        print('need arguments: path, n, save_every, timesteps')
    path = sys.argv[1]
    if path[-1] == '/':
        path = path[:-1]
    n = int(sys.argv[2])
    save_every = int(sys.argv[3])
    
    times = []
    for arg in sys.argv[4:]:
        times.append(int(arg))
    
    if len(times) < 2:
        colors = ['blue']
    elif len(times) < 3:
        colors = ['red','blue']
    else:
        colors = plot_colors()
    
    r = np.arange(n)
    
    for i in range(len(times)):
        t = times[i]
        psi = read_wave(path, t, save_every)
        rho = (np.abs(psi) * r)**2
        plt.plot(r, rho, '.', color = colors[i % len(colors)],
                    label = 'step %d' % t)
    plt.legend()
    plt.show()