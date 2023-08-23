# -*- coding: utf-8 -*-

# This code is part of Schroedinger Newton Inhibitions of Dispersion
# (C) Copyright Andre Grossardt 2010-2023
# https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
# 
# This code is licensed under the MIT License (see LICENSE.txt for details)

"""
This script is used to create a movie of phase of the wave function
created by the simulation.

Called with the following arguments (same syntax as movie.py):

phase.py path n max_t save_every [plot_title w m dr dt]

    path: path to the data files
    n: number of grid points
    max_t: maximum time up to which the data is read
    save_every: save every nth time step (parameter of simulation)
    plot_title: title of the plot (optional)
    w: width of the initial gaussian (optional)
    m: mass of the particle (required if w is given)
    dr: grid step size (required if w is given)
    dt: time step size (required if w is given)

alternative call (Plots two solutions in one plot):
movie.py path n max_t save_every plot_title second_path
    
    second_path: path to the data files of the second solution
"""

import sys
import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from helpers import free_solution, read_waves, phase


if __name__ == "__main__":
    # arguments: path, n, max_t, save_every, [plot_title, w, m, dr, dt]
    # or: path, n, max_t, save_every, plot_title, second_path
    if len(sys.argv) < 5:
        print('need arguments: path, n, max_t, save_every')
        print('optional: plot_title, w, m, dr, dt')
        exit()
    path = sys.argv[1]
    if path[-1] != '/':
        path = path + '/'
    if not os.path.exists(os.path.dirname(path)):
        print('Error: Path does not exist')
        exit()
    n = int(sys.argv[2])
    max_t = int(sys.argv[3])
    save_every = int(sys.argv[4])
    out_path = path + 'movie_tmp/'
    movie_file = path + 'phase_mov.avi'
    if not os.path.exists(os.path.dirname(out_path)):
        os.makedirs(os.path.dirname(out_path))
    
    if len(sys.argv) > 5:
        plot_title = sys.argv[5]
    else:
        plot_title = ""
    plot_free = False
    plot_second = False
    if len(sys.argv) > 9:
        plot_free = True
        w = float(sys.argv[6])
        m = float(sys.argv[7])
        dr = float(sys.argv[8])
        dt = float(sys.argv[9])
    elif len(sys.argv) > 6:
        plot_second = True
        second_path = sys.argv[6]
        if second_path[-1] != '/':
            second_path = second_path + '/'
        if not os.path.exists(os.path.dirname(second_path)):
            print('Error: Path does not exist')
            exit()
    
    # frames per second:
    fps = 25
    
    psi = read_waves(path, n, max_t, save_every)
    if plot_second:
        psi_sec = read_waves(second_path, n, max_t, save_every)

    # Create video, this code is mostly from
    # http://matplotlib.sourceforge.net/examples/animation/movie_demo.html
    
    not_found_msg = """
    The mencoder command was not found;
    mencoder is used by this script to make an avi file from a set of pngs.
    It is typically not installed by default on linux distros because of
    legal restrictions, but it is widely available.
    """

    try:
        subprocess.check_call(['mencoder'])
    except subprocess.CalledProcessError:
        print('mencoder command was found')
        pass # mencoder is found, but returns non-zero exit as expected
        # This is a quick and dirty check; it leaves some spurious output
        # for the user to puzzle over.
    except OSError:
        print(not_found_msg)
        sys.exit("quitting\n")

    if plot_free:
        x = np.arange(n) * dr
    else:
        x = np.arange(n)

    for i in range(len(psi)) :
        if plot_free:
            yf = phase(free_solution(w, m, (i+1) * save_every * dt, n, dr))
            plt.plot(x,yf,'r.', label='free')
        elif plot_second:
            plt.plot(x,yf[i],'r.', label='free')
        plt.plot(x,phase(psi[i]),'b.', label='grav.')
        plt.axis((x[0],x[-1],-3.2,3.2))
        plt.xlabel('r (m)')
        plt.ylabel('phase')
        plt.legend()
        plt.title(plot_title, fontsize=20)

        filename = out_path + str('%08d' % i) + '.png'
        plt.savefig(filename, dpi=100)
        print('Wrote file', filename)
        # Clear the figure to make way for the next image.
        plt.clf()

    # emulate console call of mencoder
    # mencoder mf://*.png -mf type=png:w=800:h=600:fps=25 -ovc lavc -lavcopts
    #          vcodec=mpeg4 -oac copy -o output.avi
    # See the MPlayer and Mencoder documentation for details.

    command = ('mencoder',
               'mf://' + out_path + '*.png',
               '-mf',
               'type=png:w=800:h=600:fps=' + str(fps),
               '-ovc',
               'lavc',
               '-lavcopts',
               'vcodec=mpeg4',
               '-oac',
               'copy',
               '-o',
               movie_file)

    print("\n\nabout to execute:\n%s\n\n" % ' '.join(command))
    subprocess.check_call(command)

    print("\n\n Deleting %s\n\n" % out_path)
    shutil.rmtree(out_path)

    print("\n\n The movie was written to '%s'" % movie_file)
