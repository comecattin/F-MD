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
      \right] +
      C q_i q_j \cfrac{e^{-\alpha R_{ij}}}{R_{ij}}
  \right\}
$$

Parameters are given in the following table:
|$k_b$ (kcal mol⁻¹ Å⁻²)|$r_\text{OH}^\text{eq}$ (Å) | $k_a$ (kcal mol⁻¹ rad⁻²) | $\theta_\text{HOH}^\text{eq}$ (deg) | $\sigma_\text{OO}$ (Å) | $\epsilon_\text{OO}$ (kcal.mol-1) | $q(\text{O})$ (e) | $q(\text{H})$ (e) |
|---------|---------|--------|---------|--------|--------|--------|--------|
| 1059.162 | 1.012 | 75.90 | 113.24 | 3.165492 | 0.1554253 | -0.82 | 0.41 |

- Bonded and angle interactions are modeled by harmonic potentials
- van der Waals interactions are modeled by a Lennard-Jones potential
- Coulomb interactions are modeled by a Yukawa potential, *ie.* screened Coulomb potential

The Yukawa energy uses the positive prefactor `q_i*q_j`, so like charges repel
and the forces equal the negative energy gradient. Electrostatic pairs currently
include atoms within the same molecule as well as atoms in different molecules.

## Units

All inputs and calculations use the following units:

| Quantity | Unit |
|---|---|
| Positions, box length, bond lengths | Å |
| Time and `dt` | ps (1 ps = 1000 fs) |
| Velocities | Å/ps |
| Masses | g/mol |
| Charges | elementary charge, e |
| Potential, kinetic and total energies | kcal/mol (sum over the simulated system) |
| Forces | kcal mol⁻¹ Å⁻¹ |
| Angles in energy/force calculations | radians |
| Yukawa screening coefficient `alpha` | Å⁻¹ (currently 1.0) |

The shared constants in `src/units.f90` convert force divided by mass to
acceleration in Å/ps² with a factor of **418.4**. Kinetic energy is therefore
`sum(m*v²)/(2*418.4)` in kcal/mol. Electrostatic energy is
`C*q_i*q_j*exp(-alpha*r)/r`, with **C ≈ 332.063713 kcal Å mol⁻¹ e⁻²**,
derived from [NIST CODATA 2022 constants](https://physics.nist.gov/cuu/pdf/all.pdf).
The screening multiplier in the Yukawa helper routines is dimensionless; it is
not an atomic mass and is set to one in the simulation.

These conversions change the dynamics and energy values from earlier versions,
which omitted them. Existing timestep values are now interpreted as ps. The
default timestep is `0.0001` ps (0.1 fs); the example uses `0.00001` ps (0.01 fs).
Check timestep convergence for each simulation. Initial velocity components are
uniformly distributed between -0.5 and 0.5 Å/ps; they do not specify a temperature.

## Installation

To install, clone this repository:

```bash
git clone git@github.com:comecattin/F-MD.git
cd F-MD
```

Create and activate the Python environment with [uv](https://docs.astral.sh/uv/):

```bash
uv sync --locked
source .venv/bin/activate
```

Dependencies are declared in `pyproject.toml`; `uv.lock` records the resolved
versions. The environment includes NumPy and Matplotlib for plotting and Pillow
for GIF export. Install the `ffmpeg` executable separately for MP4 export (for
example, `sudo apt install ffmpeg` on Debian/Ubuntu).
Keep the environment active when running the simulator, which calls the Python
input parser, and when running the plotting scripts. Run `uv sync` again after
dependency changes. You can also use `uv run python ...` without activating it.

On Debian/Ubuntu, install native prerequisites with
`sudo apt install gfortran make ffmpeg`. Compile and install from the repository
root (no administrator privileges needed for this prefix):

```bash
make -j2
make install PREFIX="$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
```

This installs `md_simulation`, `fmd-plot-energies`, and `fmd-plot-trajectories` in
`PREFIX/bin`, with the executable and Python helpers in `PREFIX/libexec/f-md`.
Keep the uv environment active, or set `FMD_PYTHON` to the absolute path of its
Python executable. `DESTDIR` is supported for staged installations. Source-tree
runs with `cd src` and `./md_simulation input.in` also remain supported.

## Run an MD simulation

Run the installed simulator in the directory where outputs should be written:

```bash
mkdir -p run
cd run
md_simulation ../examples/example_input.in
```

Input paths may contain spaces. Invalid atom counts, strides, timesteps, and box
lengths fail before the simulation starts.

The file `input.in` is the input file and contain:

- The number of particles (`n_atoms`, default `30`)
- The number of time steps (`n_steps`, default `1000`)
- The time step (`dt`, in ps, default `0.0001`)
- The box length (`box_length`, in Å, default `10.0`)
- Stride, write position, energies and print them every nth step (`stride`, default `1`)

An example is given in the `examples/example_input.in` file.

## Output and plot the trajectories

By default, positions along time steps are outputted in the `trajectories.dat` file.

To plot an animation of the trajectories, run `fmd-plot-trajectories` in the directory where the `trajectories.dat` file is written.

## Output and plot the energies

Kinetic, potential and total energies are by default computed at each time step and saved under the `energies.dat` file.

To plot theses energies along the time, run `fmd-plot-energies energies.dat`. An `energies.pdf` file will be written.

On a machine without a display, use
`MPLBACKEND=Agg fmd-plot-energies energies.dat` to save the PDF.

## Physics regression checks

Run `bash tests/run.sh` from the repository root (requires `gfortran`). The checks
compare forces with numerical energy gradients, test periodic translation and
zero net internal force, and verify energy conservation and timestep convergence
for a deterministic two-water simulation. Independent physical unit checks cover
kinetic energy, electrostatic energies and forces, and both integration half-steps.
Output checks cover energies too large for fixed-width decimal formatting.
Build files are kept in a temporary directory.


## Installation and code tests

From the repository root, run the same test command used by CI:

```bash
uv sync --locked
uv run --locked make test
```

This runs the Fortran regression checks plus Python tests for configuration
validation, installation into a temporary prefix, simulations from a different
working directory, finite output, and plotting. Tests use temporary directories;
they do not overwrite simulation results. The installed MP4 export check skips
when FFmpeg is absent locally. To require it as CI does, use:

```bash
FMD_REQUIRE_FFMPEG=1 uv run --locked make test
```

GitHub Actions runs on pushes and pull requests with both optimized and
runtime-checked Fortran builds. It installs native dependencies, syncs `uv.lock`,
compiles in parallel, and runs the complete suite including GIF/MP4 export.
