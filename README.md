# F-MD, Fortran Molecular Dynamic
A simple Molecular Dynamic implementation in Fortran.

## Installation
To install, clone this repository:

```
git clone git@github.com:comecattin/F-MD.git
cd F-MD
```

Then compile the source:
```
cd src
make
```
## Run an MD simulation
To run an MD simulation simply run `./md_simulation` in the `src` directory.

Parameters for the MD simulation can be changed inside the `src/main.f90` file. Compilation is needed after making changes.

## Ouput and plot the trajectories
To plot an animation of the trajectories, run `python plot_trajectories.py` in the `src` directory. 
