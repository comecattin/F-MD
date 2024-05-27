# F-MD, Fortran Molecular Dynamic
<p align='center'>
  <img src=https://github.com/comecattin/F-MD/assets/75748278/1d46d2ed-7f07-4c99-9086-a80ddbaf25b5>
</p>

A simple Molecular Dynamic implementation in Fortran.

Molecules are randomly initialized in a simulation box with Periodic Bondary Condition using minimal image. Lennard-Jones force are then propagated along the time using a Velocity Verlet integrator.

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

An example is given in the `examples/example_input.dat` file.

## Output and plot the trajectories

By default, positions along time steps are outputted in the `trajectories.dat` file.

To plot an animation of the trajectories, run `python plot_trajectories.py` in the directory where the `trajectories.dat` file is written.
