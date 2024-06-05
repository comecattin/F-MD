#!/usr/bin/env python3
"""Plot the energies of the system as a function of time."""
import matplotlib.pyplot as plt
import numpy as np

def load_data(filename):
    """Load the data from the file."""
    with open(filename, 'r') as f:
        data = f.readlines()
    for i, line in enumerate(data):
        data[i] = float(line.replace('Step:', '').replace('Kinetic Energy:', '').replace('Potential Energy:', '').replace('Total Energy:', '').replace('\n', '').replace(' ', ''))
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

if __name__ == "__main__":
    filename = '/home/ccattin/dev/fortran_MD/src/energies.dat'
    kinetic_energy, potential_energy, total_energy, step = load_data(filename)
    plot_energies(kinetic_energy, potential_energy, total_energy, step)