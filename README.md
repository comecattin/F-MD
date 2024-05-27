# F-MD, Fortran Molecular Dynamic

A simple Molecular Dynamic implementation in Fortran.

## Installation

To install, clone this repository:

```bash
git clone git@github.com:comecattin/F-MD.git
cd F-MD
```

Then compile the source:

```bash
cd src
make
```

## Run an MD simulation

To run an MD simulation simply run `./md_simulation input.dat` in the `src` directory.

The file `input.dat` is the input file and contain in the following order:

- The number of particles
- The number of time steps
- The time step
- The box length

An example is given in the `example_input.dat` file.

## Output and plot the trajectories

By default, positions along time steps are outputted in the `trajectories.dat` file.

To plot an animation of the trajectories, run `python plot_trajectories.py` in the `src` directory.
