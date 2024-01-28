# Schrödinger-Newton Inhibitions of Dispersion

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
](https://opensource.org/licenses/MIT)

**This repository is for showcasing only. You are free to use any of its
contents, but it is currently not maintained and not open for contributions.
This code was originally not intended for publication, is poorly documented,
and may not run out of the box on your device (my apologies).**

Simulation of the dynamics of the spherically symmetric single particle
Schrödinger-Newton equation. The code uses a Crank-Nicolson scheme for the
evolution of the time dependent nonlinear Schrödinger-Newton equation and is
inspired by the PhD thesis by Salzman[^salzman2005]. It has been used to
retrieve the results published as part of my PhD research[^giulini2011].

[^salzman2005]: Salzman, Peter Jay. 2005. Investigation of the Time Dependent
Schrödinger-Newton Equation. PhD Thesis, Davis: University of California.

[^giulini2011]: Giulini, Domenico, and André Großardt. 2011. Gravitationally
Induced Inhibitions of Dispersion According to the Schrödinger-Newton Equation.
[Classical and Quantum Gravity 28 (19):
195026](https://doi.org/10.1088/0264-9381/28/19/195026).
Preprint at [arXiv:1105.1921](https://arxiv.org/abs/1105.1921).

## Installation

Requirements:

* GNU C compiler installed and gcc executable in PATH

Before running:

* Edit the path in line 22 of `run.sh` to reflect your preferred output
  directory

## Usage

The simulation is run using the `run.sh` script. The parameters for the
simulation are set in the `param.h` header file. Options for running using the
`run.sh` script provided are:

* `./run.sh s`                starts a new calculation
* `./run.sh c path time`      continues from given path and timestep
* `./run.sh safe s`           compile in safe mode with fewer optimization
* `./run.sh safe c path time` same but continues calculation

For the evaluation of the results, the following scripts are provided
in the `python_scripts` directory:

* `movie.py path n max_t save_every [plot_title w m dr dt]`
  creates a movie of the probability density
* `phase.py path n max_t save_every [plot_title w m dr dt]`
  creates a movie of the wave function phase
* `plot.py path n save_every times`
  creates a plot of the probability density for given time steps
* `plotfile.py path outfile-prefix n save_every dr times`
  writes the probability density to a file for a given time steps
  for further processing
* `plotphasefile.py path outfile-prefix n save_every dr times`
  writes the phase to a file
* `rmax.py path n max_t save_every w m dr dt outfile`
  writes a file containing the maximum radius of the wave function for each
  time step
* `r90.py path n max_t save_every w m dr dt outfile`
  writes a file containing the radius at which 90% of the probability density
  is contained
* `masstime.py path n max_t save_every w m dr dt deviation`
  prints the time at which a certain relative deviation in width of the
  probability distribution from the free evolution is reached

The parameters for these Python scripts are:

* `path`: path to the directory containing the simulation results
* `n`: number of grid points
* `max_t`: maximum time step
* `save_every`: save every nth time step (parameter from simulation run)
* `plot_title`: title of the plot
* `w`: width of the initial Gaussian wave packet
* `m`: mass of particle
* `dr`: spatial grid size
* `dt`: temporal step size
* `times`: time steps for which to plot, separated by spaces
* `outfile`: file to write to
* `outfile-prefix`: prefix for output file name (time step is appended)
* `deviation`: relative deviation in width of probability distribution from
  free evolution

## Background

We solve the the time dependent Schrödinger-Newton equation

```math
i\hbar \frac{\partial}{\partial t}\Psi(t,\vec x)=
H[\Psi] \Psi(t,\vec x).
```

Due to the nonlinearity of the equation, the Hamilton operator $H$ has a
functional dependence on $\Psi$ and is given by

```math
H[\Psi]=
-\frac{\hbar^2}{2m}\Delta
+m \Phi[\Psi]
=
-\frac{\hbar^2}{2m}\Delta
-Gm^2\int\frac{\vert\Psi(t,\vec y)\vert^2}{\Vert\vec x-\vec y\Vert}\,d^3y
```

where $\Delta$ is the Laplace operator, $G$ is the gravitational constant, $m$
the mass of the particle and $\Psi(t,\vec x)$ the wave function.
$\vert\cdot\vert$ denotes the complex modulus and $\Vert\cdot\Vert$ the
Euclidean norm in $\mathbb{R}^3$.

We consider the spherically symmetric case, i.e. $\Psi(t,\vec x)=\Psi(t,r)$,
in which case the potential $\Phi$ is given by

```math
\Phi[\Psi](r) = -4 \pi G m \left( \frac{1}{r} \int_0^r \vert\Psi(t,r')\vert^2 r'^2 \, dr'
     + \int_r^\infty \vert\Psi(t,r')\vert^2 r' \, dr' \right).
```

We define a spatial and temporal grid size $\Delta r$ and $\Delta t$ and use
the index notation $\Psi^n_j = \Psi(n \Delta t, j \Delta r)$. We can write the
Schrödinger-Newton equation in a discretized way using Cayley's form

```math
\exp \left(\frac{i \, \Delta t}{2 \hbar} H[\Psi]\right) \, \Psi^{n+1}_j
= \exp \left(-\frac{i \, \Delta t}{2 \hbar} H[\Psi]\right) \, \Psi^{n}_j.
```

Linearizing this equation we can write it as

```math
\Psi^{n+1}_j = (Q^{-1} - I) \Psi^{n}_j
```

with the matrix

```math
Q = \frac{1}{2} \, \left(I
     + \frac{i \, \Delta t}{2 \hbar} \, H\right),
```

where $I$ is the identity matrix. We thus have to solve the linear system

```math
Q \, \chi^n = \Psi^n
```

in order to obtain the wave function at the next time step:

```math
\Psi^{n+1} = \chi^n - \Psi^n.
```

The radial component of the Laplacian in spherical coordinates is

```math
\Delta_r = \left\{\begin{array}{ll} \frac{\partial^2}{\partial r^2}
           + \frac{2}{r} \frac{\partial}{\partial r}
           & \quad \text{if } r > 0\\ 3 \frac{\partial^2}{\partial r^2}
           &  \quad \text{if } r = 0, \end{array}\right.
```

It takes the discretized form

```math
\Delta \chi^n_j = \left\{\begin{array}{ll} \frac{1}{(\Delta r)^2} \,
                  \left( \frac{j-1}{j}\, \chi^n_{j-1} - 2 \chi^n_{j}
                  + \frac{j+1}{j}\, \chi^n_{j+1} \right)
                  & \quad \text{if } j > 0\\ \frac{1}{(\Delta r)^2} \,
                  \left( -6 \chi^n_0 + 6 \chi^n_1 \right)
                  &  \quad \text{if } j = 0, \end{array}\right.
```

while the discretized form of the potential $\Phi$ is
$\Phi^n_j = -4 \pi G m (\Delta r)^2 v^n_j$ with

```math
v^n_j = \frac{1}{j} \sum_{i=0}^{j-1} \vert\Psi^n_i\vert^2 \, i^2
         + \sum_{i=j}^{N-1} \vert\Psi^n_i\vert^2 \, i.
```

$Q$ then becomes a tridiagonal matrix

```math
Q = \begin{pmatrix}
b_0 & c_0 & 0   & 0 & \cdots \\
a_1 & b_1 & c_1 & 0 & \cdots \\
0   & a_2 & b_2 & c_2 &  \\
\vdots &  &     & \ddots  & \vdots \\
0   & \cdots &  & a_{N-1} & b_{N-1}
\end{pmatrix}.
```

Using the shorthand notations

```math
\beta = -\frac{i \hbar}{8 m} \, \frac{\Delta t}{(\Delta r)^2} \qquad,
\gamma = \frac{i \pi G}{\hbar} \, m^2 \, \Delta t \, (\Delta r)^2,
```

the diagonal elements are given by

```math
a_j = \beta \, \frac{j-1}{j} \quad (0 < j \leq N-1),
```

and the off-diagonal elements by

```math
b_0 = \frac{1}{2} - 6 \beta - \gamma\, v_0 \,,\qquad
b_j = \frac{1}{2} - 2 \beta - \gamma\, v_j  \quad (0 < j \leq N-1)
```

and

```math
c_0 = 6 \beta \,, \qquad
c_j = \beta \, \frac{j+1}{j} \quad (0 < j < N-1)
```

We use the tridiagonal matrix algorithm to solve the linear system.

## Description of files

Routines for the simulation, written in C, and bash scripts for compiling
and running the simulation are located in the `code` directory:

* `run.sh`: Bash script to compile and run the simulation
* `compile.sh`: Bash script to only compile without running
* `param.h`: Header file containing all parameters to be set
* `sne.c`: Main program running Crank-Nicolson algorithm
* `wf.c`: Definition of different wave function shapes
* `helpers.c`: Various helper functions for handling of files and output

Python scripts for evaluation of the results are located in the
`python_scripts` directory, with descriptions given above in the section
"Usage". The file `helpers.py` contains helper functions for the Python
scripts.

## License

[MIT License](LICENSE.txt)
