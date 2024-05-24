import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.animation import FuncAnimation, PillowWriter

def read_trajectories(filename):
    steps = []
    positions = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        step = None
        pos = []
        for line in lines:
            if line.startswith('Step'):
                if step is not None:
                    steps.append(step)
                    positions.append(np.array(pos))
                step = int(line.split()[1])
                pos = []
            else:
                pos.append([float(x) for x in line.split()])
        if step is not None:
            steps.append(step)
            positions.append(np.array(pos))
    return steps, positions

def animate(i):
    ax.clear()
    pos = positions[i]
    ax.scatter(pos[:, 0], pos[:, 1], s=10)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title(f'Time Step: {steps[i]}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

def main():
    global steps, positions, ax
    steps, positions = read_trajectories('trajectories.dat')

    steps = steps[::20]
    positions = positions[::20]

    fig, ax = plt.subplots(figsize=(8, 8))
    ani = FuncAnimation(fig, animate, frames=len(steps), repeat=False)

    # Save the animation
    writer = PillowWriter(fps=24)
    ani.save('particle_trajectories.gif', writer=writer)

if __name__ == '__main__':
    main()
