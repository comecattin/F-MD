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

To run an MD simulation simply run `./md_simulation input.in` in the `src` directory.

The file `input.in` is the input file and contain:

- The number of particles (`n_atoms`, default `30`)
- The number of time steps (`n_steps`, default `1000`)
- The time step (`dt`, default `0.001`)
- The box length (`box_length`, default `10.0`)
- The tolerance of the SHAKE algorithm for angle and bonds length constraints (`tolerance`, default `1e-6`)
- The maximum number of iteration for the SHAKE algorithm (`max_iter`, default `100`)

An example is given in the `examples/example_input.in` file.

## Output and plot the trajectories

By default, positions along time steps are outputted in the `trajectories.dat` file.

To plot an animation of the trajectories, run `python plot_trajectories.py` in the directory where the `trajectories.dat` file is written.

## Output and plot the energies

Kinetic, potential and total energies are by default computed at each time step and saved under the `energies.dat` file.

To plot theses energies along the time, run `python plot_energies.py energies.dat`. An `energies.pdf` file will be written.
