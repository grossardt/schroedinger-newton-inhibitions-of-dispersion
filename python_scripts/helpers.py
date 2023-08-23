# -*- coding: utf-8 -*-

# This code is part of Schroedinger Newton Inhibitions of Dispersion
# (C) Copyright Andre Grossardt 2010-2023
# https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
# 
# This code is licensed under the MIT License (see LICENSE.txt for details)

"""
Helper functions for python scripts used to create movies and plots
"""

import numpy as np
import cmath
from math import pi, sqrt

def free_solution(w, m, t, n, dr):
    """
    Exactly solve the free particle

    Args:
        w: width of the initial gaussian
        m: mass of the particle
        t: time
        n: number of points
        dr: grid step size

    Returns:
        wave function at time t
    """
    # index  [:n] is neccessary to make sure that psi has the right length
    r = np.arange(0, n * dr, dr, 'complex')[:n]
    # z = m / (m + i hbar alpha t) (dimensionless)
    z = 1. / (1. + 63.50779875974j / m / w**2 * t)
    psi = (z / sqrt(pi) / w)**1.5 * np.exp(-r**2 * z / 2 / w**2)
    return psi

def read_waves(path, n, max_t, save_every):
    """
    Read the wave function data from the files

    Args:
        path: path to the data files
        n: number of grid points
        max_t: maximum time up to which the data is read
        save_every: save every nth time step (parameter of simulation)

    Returns:
        array of wave functions
    """
    saved = max_t / save_every
    r = range(saved)
    t = (np.array(r) + 1) * save_every
    if max_t % save_every > 0:
        r.append(saved)
        t.append(max_t)
        saved += 1
    res = np.zeros((saved, n), 'complex256')
    path += "data/"
    for i in r:
        filename = '%sw%014d.dat' % (path, t[i])
        res[i] = np.fromfile(filename, 'complex256')
        if not np.isfinite(res[i]).all():
            print('ERROR: NaN in img %s' % i)
    return res

def read_wave(path, t, save_every = 1):
    """
    Reads the wave function at a single time step

    Args:
        path: path to the data files
        t: time step
        save_every: save every nth time step (parameter of simulation)
                    default: 1

    Returns:
        wave function at time t
    """
    filename = '%s/data/w%014d.dat' % (path, t * save_every)
    return np.fromfile(filename, 'complex256')

def phase(psi):
    """
    Returns the phase of a given wave function

    Args:
        psi: wave function

    Returns:
        phase of the wave function
    """
    p = np.zeros(len(psi))
    for i in range(len(psi)):
        p[i] = cmath.phase(psi[i])
    return p

def html_color(r, g, b):
    """
    Returns a html color string from rgb values

    Args:
        r: red value
        g: green value
        b: blue value

    Returns:
        html color string
    """
    string = '#'
    for x in r, g, b:
        if x < 16:
            string += '0'
        string += hex(x)[2:]
    return string

def plot_colors():
    """
    Returns a list of html colors for plotting
    """
    return [html_color(0,72,140), html_color(32,95,154),
            html_color(62,117,167), html_color(119,157,191),
            html_color(178,200,217), html_color(181,181,169),
            html_color(138,141,130),  html_color(116,121,111),
            html_color(93,100,90), html_color(67,76,67),
            html_color(51,38,0), html_color(99,82,43),
            html_color(122,101,62), html_color(144,122,82),
            html_color(190,163,122), html_color(169,204,143),
            html_color(130,177,106), html_color(92,151,70),
            html_color(61,129,40), html_color(30,108,11),
            html_color(240,180,0), html_color(243,192,28),
            html_color(248,215,83), html_color(250,225,107),
            html_color(255,248,163), html_color(179,0,35),
            html_color(188,28,57), html_color(196,56,79),
            html_color(214,112,123), html_color(230,165,164)
            ]