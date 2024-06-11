# F-MD, Fortran Molecular Dynamic

<p align='center'>
  <img src=https://github.com/comecattin/F-MD/assets/75748278/1d46d2ed-7f07-4c99-9086-a80ddbaf25b5>
</p>

A simple Molecular Dynamic implementation in Fortran.

Water molecules are randomly initialized in a simulation box with Periodic Bondary Condition using minimal image. This work relies on the SPC/Fw water model, where interaction potential is described as:

$$V = V^\text{intra} + V^\text{inter}$$

with
$$
  V^\text{intra} = \cfrac{k_b}{2}\left[
    (r_{\text{OH}_1} - r_\text{OH}^\text{eq})^2
    + (r_{\text{OH}_2} - r_\text{OH}^\text{eq})^2
    \right]
    + \cfrac{k_a}{2} (\theta_\text{HOH} - \theta_\text{HOH}^\text{eq})^2
$$

and
$$
  V^\text{inter} = \sum_{ij}^\text{all pairs} \left\{
    4\epsilon_{ij} \left[
      \left(\cfrac{\sigma_{ij}}{R_{ij}} \right)^{12} -
      \left(\cfrac{\sigma_{ij}}{R_{ij}} \right)^{6}
      \right] -
      q_i q_j \cfrac{e^{-\alpha m R_{ij}}}{R_{ij}}
  \right\}
$$

Parameters are given in the following table:
|$k_b$|$r_\text{OH}^\text{eq}$ (Å) | $k_a$ | $\theta_\text{HOH}^\text{eq}$ (deg) | $\sigma_\text{OO}$ (Å) | $\epsilon_\text{OO}$ (kcal.mol-1) | $q(\text{O})$ (e) | $q(\text{H})$ (e) |
|---------|---------|--------|---------|--------|--------|--------|--------|
| 1059.162 | 1.012 | 75.90 | 113.24 | 3.165492 | 0.1554253 | -0.82 | 0.41 |

- Bonded and angle interactions are modeled by harmonic potentials
- van der Waals interactions are modeled by a Lennard-Jones potential
- Coulomb interactions are modeled by a Yukawa potential, *ie.* screened Coulomb potential

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
- Stride, write position, energies and print them every nth step (`stride`, default `1`)

An example is given in the `examples/example_input.in` file.

## Output and plot the trajectories

By default, positions along time steps are outputted in the `trajectories.dat` file.

To plot an animation of the trajectories, run `python plot_trajectories.py` in the directory where the `trajectories.dat` file is written.

## Output and plot the energies

Kinetic, potential and total energies are by default computed at each time step and saved under the `energies.dat` file.

To plot theses energies along the time, run `python plot_energies.py energies.dat`. An `energies.pdf` file will be written.
