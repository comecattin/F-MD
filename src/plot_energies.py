#!/usr/bin/env python3
"""Plot the energies of the system as a function of time."""
import matplotlib.pyplot as plt
import numpy as np
import argparse

def load_data(filename):
    """Load the data from the file."""
    with open(filename, 'r') as f:
        data = f.readlines()
    for i, line in enumerate(data):
        try:
            data[i] = float(line.replace('Step:', '').replace('Kinetic Energy:', '').replace('Potential Energy:', '').replace('Total Energy:', '').replace('\n', '').replace(' ', ''))
        except ValueError:
            data[i] = None
    step = np.array(data[0::4], dtype=int)
    kinetic_energy = np.array(data[1::4])
    potential_energy = np.array(data[2::4])
    total_energy = np.array(data[3::4])
    return kinetic_energy, potential_energy, total_energy, step

def plot_energies(kinetic_energy, potential_energy, total_energy, step):
    """Plot the energies as a function of time."""
    plt.plot(step, kinetic_energy, label='Kinetic Energy')
    plt.plot(step, potential_energy, label='Potential Energy')
    plt.plot(step, total_energy, label='Total Energy')
    plt.xlabel('Time step')
    plt.ylabel('Energy')
    plt.legend()
    plt.savefig('energies.pdf', bbox_inches='tight', dpi=300)
    plt.show()

def main():
    """Read the filename from the command line and plot the energies."""
    parser = argparse.ArgumentParser(
        description='Plot the energies of the system as a function of time.'
    )
    parser.add_argument(
        'filename',
        help='The name of the file containing the energies.'
    )
    args = parser.parse_args()
    kinetic_energy, potential_energy, total_energy, step = load_data(args.filename)
    plot_energies(kinetic_energy, potential_energy, total_energy, step)

if __name__ == "__main__":
    main()